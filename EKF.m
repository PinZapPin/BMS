clc;
clearvars;
close all;

%% 1. MEMBACA DATA DARI FILE CSV
try
    rc_params_time_varying = readtable('Estimated_ECM_Parameters_and_SOC.csv');
    measurement_data = readtable('DSTEKF.csv'); 
catch ME
    error('Pastikan file CSV ada di direktori yang sama. Error: %s', ME.message);
end

%% 1.1. PARAMETER STATIS
eta = 0.95; 
C_nominal = 1.085; % Kapasitas nominal dalam Ah

Q_val_v1 = 4e-6;  
Q_val_v2 = 4e-6;  
Q_val_s = 4e-5;   
R_meas = 0.02;

ocv_coeffs_cha = [3813; -16906; 31678; -32618; 20076; -7522; 1671; -205; 12; 3];
ocv_coeffs_dis = [3279; -15807; 32591; -37545; 26475; -11764; 3265; -541; 49; 1];

alpha_ocv = 20; % Parameter transisi histeresis OCV

%% 2. PERSIAPAN DATA PENGUKURAN
time_data    = measurement_data.Time_s;
current_data = measurement_data.Current_A;
voltage_data = measurement_data.Voltage_V;

dt = time_data(2) - time_data(1);
N  = length(time_data);

R0_data = rc_params_time_varying.R0;
R1_data = rc_params_time_varying.R1;
R2_data = rc_params_time_varying.R2;
C1_data = rc_params_time_varying.C1;
C2_data = rc_params_time_varying.C2;

%% 2.1. HITUNG SOC REFERENSI
true_soc_init = 100; 
true_soc = calculate_soc_coulomb_counting(current_data, time_data, C_nominal, true_soc_init);

%% 3. INISIALISASI EKF
x_est  = zeros(3, N); 
P_est  = zeros(3, 3, N);

x_est(:,1)    = [0.01; 0.01; true_soc_init/100];
P_est(:,:,1) = diag([0.1, 0.1, 0.01]);

Q = diag([Q_val_v1, Q_val_v2, Q_val_s]);
R_k = R_meas;

%% 4. EKF LOOP UTAMA
for k = 2:N
    current_k   = current_data(k);
    voltage_k   = voltage_data(k);
    current_km1 = current_data(k-1);

    x_prev = x_est(:,k-1);
    P_prev = P_est(:,:,k-1);

    R0 = max(R0_data(k), 1e-9);
    R1 = max(R1_data(k), 1e-9);
    R2 = max(R2_data(k), 1e-9);
    C1 = max(C1_data(k), 1e-9);
    C2 = max(C2_data(k), 1e-9);

    v1_new = (1 - dt/(R1*C1)) * x_prev(1) + (dt/C1) * current_km1;
    v2_new = (1 - dt/(R2*C2)) * x_prev(2) + (dt/C2) * current_km1;
    delta_soc = (eta * current_km1 * dt) / (C_nominal * 3600);
    soc_new = max(0, min(1, x_prev(3) - delta_soc));
    x_priori = [v1_new; v2_new; soc_new];

    A = [1 - dt/(R1*C1), 0, 0; 0, 1 - dt/(R2*C2), 0; 0, 0, 1];
    P_priori = A * P_prev * A' + Q;

    s = x_priori(3);

    % OCV dari masing-masing kurva
    ocv_cha = polyval(ocv_coeffs_cha, s);
    ocv_dis = polyval(ocv_coeffs_dis, s);

    % Smooth sigmoid-based blending weight
    sigma = 1 / (1 + exp(-alpha_ocv * current_k)); % sigmoid(current)

    % OCV halus (histeresis)
    ocv = sigma * ocv_cha + (1 - sigma) * ocv_dis;

    % dOCV/dSOC untuk masing-masing
    dOCVds_cha = polyder(ocv_coeffs_cha);
    dOCVds_dis = polyder(ocv_coeffs_dis);

    dOCVds_val_cha = polyval(dOCVds_cha, s);
    dOCVds_val_dis = polyval(dOCVds_dis, s);

    % dOCV/dSOC hasil blending
    dOCVds = sigma * dOCVds_val_cha + (1 - sigma) * dOCVds_val_dis;

    Ck = [-1, -1, dOCVds];
    L = P_priori * Ck' / (Ck * P_priori * Ck' + R_k);

    v_pred = ocv - x_priori(1) - x_priori(2) - current_k * R0;
    x_est(:,k) = x_priori + L * (voltage_k - v_pred);
    x_est(3,k) = max(0, min(1, x_est(3,k)));
    P_est(:,:,k) = (eye(3) - L * Ck) * P_priori;
end

%% Konversi SoC EKF ke persen
estimated_soc_percent = x_est(3,:) * 100;

%% 5. GRAFIK HASIL EKF
figure('Name', 'Hasil Final EKF', 'Position', [100, 100, 900, 500]);
plot(time_data, true_soc, 'k', 'LineWidth', 2); hold on;
plot(time_data, estimated_soc_percent, 'r--', 'LineWidth', 1.5); hold off;
legend('True SOC (Referensi)', 'Estimated SOC (EKF)');
xlabel('Time (s)'); ylabel('SOC (%)');
title('True SOC vs Estimated SOC (EKF)');
grid on; ylim([-5 105]); 
set(gca, 'FontSize', 12);

%% 6. ANALISIS ERROR
error = true_soc - estimated_soc_percent';
rmse = sqrt(mean(error.^2));
mae  = mean(abs(error));
max_error = max(abs(error));

fprintf('Analisis Error Estimasi SOC:\n');
fprintf('=============================\n');
fprintf('Root Mean Square Error (RMSE): %.4f %%\n', rmse);
fprintf('Mean Absolute Error (MAE)    : %.4f %%\n', mae);
fprintf('Maximum Absolute Error       : %.4f %%\n', max_error);

%% 7. PLOT OCV DAN TURUNANNYA
s = linspace(0, 1, 1000);
ocv_cha = polyval(ocv_coeffs_cha, s);
ocv_dis = polyval(ocv_coeffs_dis, s);

figure('Name', 'OCV vs SOC', 'Position', [200, 200, 900, 400]);
plot(s*100, ocv_cha, 'b', 'LineWidth', 2); hold on;
plot(s*100, ocv_dis, 'r', 'LineWidth', 2); hold off;
xlabel('SOC (%)'); ylabel('OCV (V)');
legend('OCV Charging', 'OCV Discharging');
title('Kurva OCV vs SOC');
grid on;

%% FUNGSI BANTU: Hitung SOC dari Coulomb Counting
function soc = calculate_soc_coulomb_counting(current, time, Q_nominal, soc_init)
    dt_local = [0; diff(time)];
    soc = zeros(length(current), 1); soc(1) = soc_init;
    for k = 2:length(current)
        delta_soc = (current(k) * dt_local(k)) / (Q_nominal * 3600) * 100;
        soc(k) = soc(k-1) + delta_soc;
    end
    soc = max(0, min(100, soc));
end
