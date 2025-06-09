clc;
clearvars;
close all;

%% 1. Membaca Data yang Berubah Seiring Waktu dari File CSV
% Pastikan file CSV berada di folder yang sama atau tentukan path lengkapnya
try
    % Data parameter model RC yang berubah seiring waktu
    rc_params_time_varying = readtable('Estimated_ECM_Parameters_and_SOC.csv'); 
    
    % Data pengukuran real-time (OC_Current_Voltage.csv digunakan sebagai input measurement_data)
    measurement_data = readtable('DSTEKF.csv'); % Ini adalah DSTEKF.csv yang sekarang berisi data OCV

    % Error handling untuk file missing
catch ME
    error('Pastikan semua file CSV yang dibutuhkan ada di direktori yang sama: parameters_time_varyingOCV.csv, OC_Current_Voltage.csv. Error: %s', ME.message);
end

%% 1.1. Definisi Parameter Statis Langsung di Program
% Parameter statis (tidak perlu dari CSV lagi)
eta = 1; % Efisiensi Coulomb
C_nominal = 1.1; % Kapasitas Nominal (Ah)

% Kovarians noise EKF (Q)
Q_val_v1 = 0.1;
Q_val_v2 = 0.1;
Q_val_s = 0.01;

% Kovarians noise pengukuran (R)
R_meas = 1000; 

% Koefisien polinomial OCV-SOC (orde 9, didefinisikan langsung di program)
% Penting: Anda harus mendapatkan koefisien ini dari hasil polyfit orde 9 pada data OCV-SOC Anda.
% Ganti dengan nilai yang sebenarnya dari hasil polyfit!

% Koefisien untuk fungsi Charge (arus positif)
ocv_coeffs_cha = [
    -3; % a1 (koefisien untuk s^9) - GANTI DENGAN NILAI ANDA
    12;  % a2 (koefisien untuk s^8)
    -205; % a3 (koefisien untuk s^7)
    1671;  % a4 (koefisien untuk s^6)
    -7522; % a5 (koefisien untuk s^5)
    20076;  % a6 (koefisien untuk s^4)
    -32618; % a7 (koefisien untuk s^3)
    31678;  % a8 (koefisien untuk s^2)
    -16906; % a9 (koefisien untuk s^1)
    3813    % a10 (konstanta, s^0)
];

% Koefisien untuk fungsi Discharge (arus negatif)
ocv_coeffs_dis = [
    1; % a1 (koefisien untuk s^9) - GANTI DENGAN NILAI ANDA
    49;  % a2 (koefisien untuk s^8)
    -541; % a3 (koefisien untuk s^7)
    3265;  % a4 (koefisien untuk s^6)
    -11764; % a5 (koefisien untuk s^5)
    26475;  % a6 (koefisien untuk s^4)
    -37545; % a7 (koefisien untuk s^3)
    32591;  % a8 (koefisien untuk s^2)
    -15907; % a9 (koefisien untuk s^1)
    3279    % a10 (konstanta, s^0)
];


%% 2. Mempersiapkan Data Pengukuran
time_data = measurement_data.Time_s; % Dari OC_Current_Voltage.csv
current_data = measurement_data.Current_A; % i(t)
voltage_data = measurement_data.Voltage_V; % y(t)

dt = time_data(2) - time_data(1); % Interval waktu sampling
N = length(time_data); % Jumlah total langkah waktu

% Pastikan data parameter RC punya jumlah baris yang sama dengan data pengukuran
if size(rc_params_time_varying, 1) ~= N
    error('Jumlah baris pada parameters_time_varyingOCV.csv harus sama dengan OC_Current_Voltage.csv.');
end

% Ekstrak data parameter RC yang berubah seiring waktu
R0_data = rc_params_time_varying.R0;
R1_data = rc_params_time_varying.R1;
R2_data = rc_params_time_varying.R2;
C1_data = rc_params_time_varying.C1;
C2_data = rc_params_time_varying.C2;


%% 2.1. Menghitung SOC Asli (True SOC) dari Coulomb Counting
% Asumsi SOC awal untuk OC_Current_Voltage.csv diketahui (misal 100%)
true_soc_init = 100; % Sesuaikan jika SOC awal baterai Anda di pengujian ini bukan 100%
true_soc = calculate_soc_coulomb_counting(current_data, time_data, C_nominal, true_soc_init);


%% 3. Inisialisasi EKF
% State vector x = [v1; v2; s]
x_est = zeros(3, N); % Menyimpan estimasi keadaan posterior EKF
P_est = zeros(3, 3, N); % Menyimpan kovarians kesalahan posterior EKF

% Inisialisasi awal EKF
initial_soc = 0.5; % Tebakan awal SOC EKF (normalized 0-1)
x_est(:,1) = [0; 0; initial_soc]; % v1, v2 diasumsikan 0 pada awal
P_est(:,:,1) = diag([0.1, 0.1, 0.01]); % Kovarians kesalahan awal

% Kovarians noise EKF (Q)
Q = diag([Q_val_v1, Q_val_v2, Q_val_s]);
% Kovarians noise EKF (R)
R_k = R_meas;


%% 4. Loop EKF Utama
for k = 2:N
    % Ambil input pengukuran pada langkah ini
    current_k = current_data(k);         % Arus pada t(k)
    voltage_k = voltage_data(k);         % Tegangan pada t(k)
    current_k_minus_1 = current_data(k-1); % Arus pada t(k-1) untuk prediksi EKF

    % State posterior EKF dari langkah sebelumnya
    x_prev = x_est(:,k-1);
    P_prev = P_est(:,:,k-1);

    % --- Membaca Parameter R & C untuk Langkah Waktu Saat Ini ---
    R0_current = R0_data(k); 
    R1_current = R1_data(k);
    C1_current = C1_data(k);
    R2_current = R2_data(k);
    C2_current = C2_data(k);
    
    % Memastikan parameter RC valid (tidak nol atau terlalu kecil, untuk menghindari inf/NaN)
    R0_current = max(R0_current, 1e-9);
    R1_current = max(R1_current, 1e-9);
    C1_current = max(C1_current, 1e-9);
    R2_current = max(R2_current, 1e-9);
    C2_current = max(C2_current, 1e-9);

    % --- Memilih Koefisien OCV-SOC Berdasarkan Arah Arus ---
    % Konvensi: Positif = Discharge, Negatif = Charge (sesuai jurnal EKF)
    if current_k >= 0 % Asumsi: Discharge (arus positif), pakai koefisien discharge
        current_ocv_coeffs = ocv_coeffs_dis; 
    else % Asumsi: Charge (arus negatif), pakai koefisien charge
        current_ocv_coeffs = ocv_coeffs_cha;
    end


    % --- FASE 1: Time Update (Prediction Step) ---
    % 1. Prediksi Keadaan (State Estimate Time Update)
    % Model diskrit A dan B EKF, kini menggunakan R,C dari data CSV di setiap langkah k
    A_k_minus_1_ekf = [(-1/(R1_current*C1_current))*dt + 1, 0, 0;
                       0, (-1/(R2_current*C2_current))*dt + 1, 0;
                       0, 0, 1];
    B_k_minus_1_ekf = [dt/C1_current; dt/C2_current; -eta*dt];

    x_priori_ekf = A_k_minus_1_ekf * x_prev + B_k_minus_1_ekf * current_k_minus_1;

    % Pastikan SOC tetap di dalam rentang 0-1
    x_priori_ekf(3) = max(0, min(1, x_priori_ekf(3)));

    % 2. Prediksi Kovarians Kesalahan (Error Covariance Time Update)
    P_priori_ekf = A_k_minus_1_ekf * P_prev * A_k_minus_1_ekf' + Q;

    % --- FASE 2: Measurement Update (Correction Step) ---

    % 3. Hitung Matriks Jacobian Pengukuran (Ck) EKF
    soc_priori_norm_ekf = x_priori_ekf(3);
    
    % Perhitungan d_f_ds (df/ds) dari polinomial orde 9 OCV-SOC (menggunakan koefisien yang dipilih)
    % Jumlah koefisien = 10 untuk orde 9
    d_f_ds_ekf = 9*current_ocv_coeffs(1)*soc_priori_norm_ekf.^8 + ...
                 8*current_ocv_coeffs(2)*soc_priori_norm_ekf.^7 + ...
                 7*current_ocv_coeffs(3)*soc_priori_norm_ekf.^6 + ...
                 6*current_ocv_coeffs(4)*soc_priori_norm_ekf.^5 + ...
                 5*current_ocv_coeffs(5)*soc_priori_norm_ekf.^4 + ...
                 4*current_ocv_coeffs(6)*soc_priori_norm_ekf.^3 + ...
                 3*current_ocv_coeffs(7)*soc_priori_norm_ekf.^2 + ...
                 2*current_ocv_coeffs(8)*soc_priori_norm_ekf + ...
                 current_ocv_coeffs(9); % Koefisien ke-9 untuk s^1

    C_k_ekf = [1, 1, d_f_ds_ekf];

    % 4. Hitung Gain Kalman (Lk) EKF
    L_k_ekf = P_priori_ekf * C_k_ekf' / (C_k_ekf * P_priori_ekf * C_k_ekf' + R_k);

    % 5. Perbarui Keadaan (State Estimate Measurement Update) EKF
    v1_priori_ekf = x_priori_ekf(1);
    v2_priori_ekf = x_priori_ekf(2);
    
    ocv_priori_ekf = polyval(current_ocv_coeffs, soc_priori_norm_ekf); % Menggunakan koefisien yang dipilih
    
    predicted_voltage_ekf = ocv_priori_ekf + v1_priori_ekf + v2_priori_ekf + current_k * R0_current;

    x_est(:,k) = x_priori_ekf + L_k_ekf * (voltage_k - predicted_voltage_ekf);

    x_est(3,k) = max(0, min(1, x_est(3,k)));

    % 6. Perbarui Kovarians Kesalahan (Error Covariance Measurement Update) EKF
    P_est(:,:,k) = (eye(3) - L_k_ekf * C_k_ekf) * P_priori_ekf;

end


%% 5. Visualisasi Hasil
figure;
% Mengubah subplot menjadi 4 baris untuk menambahkan plot perbandingan SOC
subplot(4,1,1);
plot(time_data, x_est(3,:)*100, 'b', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Estimated SOC (%)');
title('EKF Estimated SOC');
grid on;

subplot(4,1,2);
plot(time_data, true_soc, 'k', 'LineWidth', 1.2);
hold on;
plot(time_data, x_est(3,:)*100, 'r--', 'LineWidth', 1.5);
legend('True SOC (Coulomb Counting)', 'Estimated SOC (EKF)');
xlabel('Time (s)');
ylabel('SOC (%)');
title('Comparison: True SOC vs Estimated SOC');
grid on;


subplot(4,1,3);
plot(time_data, voltage_data, 'k', 'LineWidth', 1);
hold on;
% Hitung prediksi tegangan berdasarkan estimasi SOC akhir EKF
predicted_voltage_ekf_final = zeros(1, N);
for k = 1:N
    v1_est = x_est(1,k);
    v2_est = x_est(2,k);
    soc_est_norm = x_est(3,k);
    
    % Perlu mengambil R0 yang sesuai dengan waktu 'k' juga untuk visualisasi ini
    R0_at_k = R0_data(k); 
    
    % Memilih koefisien OCV-SOC untuk visualisasi berdasarkan arah arus
    if current_data(k) >= 0 % Asumsi: Discharge (arus positif), pakai koefisien discharge
        ocv_coeffs_vis = ocv_coeffs_dis; 
    else % Asumsi: Charge (arus negatif), pakai koefisien charge
        ocv_coeffs_vis = ocv_coeffs_cha;
    end
    
    ocv_est = polyval(ocv_coeffs_vis, soc_est_norm);
    
    predicted_voltage_ekf_final(k) = ocv_est + v1_est + v2_est + current_data(k) * R0_at_k;
end
plot(time_data, predicted_voltage_ekf_final, 'r--', 'LineWidth', 1.5);
legend('Measured Voltage', 'Predicted Voltage (from EKF states)');
xlabel('Time (s)');
ylabel('Voltage (V)');
title('Measured vs Predicted Voltage');
grid on;

subplot(4,1,4);
yyaxis left;
plot(time_data, R0_data, 'b', 'DisplayName', 'R0');
hold on;
plot(time_data, R1_data, 'r', 'DisplayName', 'R1');
plot(time_data, R2_data, 'g', 'DisplayName', 'R2');
ylabel('Resistance (Ohm)');
yyaxis right;
plot(time_data, C1_data/1000, 'm--', 'DisplayName', 'C1 (kF)'); % Skala C agar terlihat
plot(time_data, C2_data, 'c--', 'DisplayName', 'C2 (F)');
ylabel('Capacitance');
xlabel('Time (s)');
title('RC Parameters from CSV (Time-Varying)');
legend('show', 'Location', 'best');
grid on;


%% Fungsi pembantu untuk Coulomb Counting
function soc = calculate_soc_coulomb_counting(current, time, Q_nominal, soc_init)
    % Fungsi ini menghitung SOC menggunakan metode Coulomb Counting.
    % current: Vektor arus (A)
    % time: Vektor waktu yang sesuai (s)
    % Q_nominal: Kapasitas nominal baterai (Ah)
    % soc_init: SOC awal dalam persen (0-100)

    dt_local = [0; diff(time)]; % Perbedaan waktu antar sampel
    soc = zeros(length(current), 1); % Inisialisasi vektor SOC
    soc(1) = soc_init; % Set SOC awal

    % Loop melalui data untuk menghitung SOC
    for k = 2:length(current)
        % Perhitungan perubahan SOC per langkah waktu (dalam persen)
        % Konvensi: Positif = Discharge, Negatif = Charge
        % Jika current(k) positif (discharge), maka SOC harus berkurang
        % Jika current(k) negatif (charge), maka SOC harus bertambah
        % Jadi, delta_soc perlu disesuaikan tandanya
        delta_soc_fraction = (current(k) * dt_local(k)) / (Q_nominal * 3600); % Perubahan dalam fraksi kapasitas
        
        % Menyesuaikan tanda berdasarkan konvensi jurnal (arus positif = discharge -> SOC turun)
        soc(k) = soc(k-1) - delta_soc_fraction * 100; % Kali 100 untuk persentase
    end

    % Memastikan SOC tetap dalam rentang 0-100%
    soc = max(0, min(100, soc)); 
end