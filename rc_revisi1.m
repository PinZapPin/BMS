clc;
clearvars;

% Mencari OCV vs SOC
data = readtable('Database/OC_Current_Voltage.csv');
Q_nominal = 1.1;
time = data.Step_Time;
current = data.Current_A;
voltage = data.Voltage_V;
dt = [0; diff(time)];

idx_cha = current > 0.0;
idx_dis = current < 0.0;

current_cha = current(idx_cha);
voltage_cha = voltage(idx_cha);
time_cha = time(idx_cha);
dt_cha = dt(idx_cha);

current_dis = current(idx_dis);
voltage_dis = voltage(idx_dis);
time_dis = time(idx_dis);
dt_dis = dt(idx_dis);

soc_cha_init = 0; 
soc_dis_init = 100;

soc_cha = calculate_soc(current_cha, time_cha, Q_nominal, soc_cha_init);   
soc_dis = calculate_soc(current_dis, time_dis, Q_nominal, soc_dis_init); 

soc_dis_norm = soc_dis / 100;
p_dis = polyfit(soc_dis_norm, voltage_dis, 9);
soc_dis_fit = linspace(min(soc_dis_norm), max(soc_dis_norm), 200);
voltage_dis_fit = polyval(p_dis, soc_dis_fit);

soc_cha_norm = soc_cha / 100;
p_cha = polyfit(soc_cha_norm, voltage_cha, 9);
soc_cha_fit = linspace(min(soc_cha_norm), max(soc_cha_norm), 200);
voltage_cha_fit = polyval(p_cha, soc_cha_fit);

% Menghitung SOC dari data DST
data_dst = readtable('Database/DST.csv');
time_dst = data_dst{:,1};

% figure;
% num_cols = width(data_dst);
% for i = 2:num_cols
%     subplot(num_cols-1,1,i-1);
%     plot(time_dst, data_dst{:,i});
%     xlabel('Time (s)');
%     ylabel(data_dst.Properties.VariableNames{i}, 'Interpreter', 'none');
%     title(['DST: ', data_dst.Properties.VariableNames{i}, ' vs Time']);
%     grid on;
% end

% Ambil data arus dan tegangan dari data DST
current_dst = data_dst.Current_A;
voltage_dst = data_dst.Voltage_V;
dt_dst = [0; diff(time_dst)];

% Hitung SOC dari coulomb counting
soc_dst_init = 100; 
soc_dst = calculate_soc(current_dst, time_dst, Q_nominal, soc_dst_init);

% Hitung Uoc dari fungsi polinomial yang sesuai arah arus
uoc_dst = zeros(size(soc_dst));
for k = 1:length(soc_dst)
    if current_dst(k) >= 0
        uoc_dst(k) = polyval(p_cha, soc_dst(k)/100); % fungsi charge
    else
        uoc_dst(k) = polyval(p_dis, soc_dst(k)/100); % fungsi discharge
    end
end

figure;

% Baris 1: Tegangan dan Arus terhadap waktu
subplot(3,1,1);
yyaxis left;
plot(time_dst, voltage_dst, 'b');
ylabel('Voltage (V)');
yyaxis right;
plot(time_dst, current_dst, 'r');
ylabel('Current (A)');
xlabel('Time (s)');
title('Voltage and Current');
grid on;

% Baris 2: SOC terhadap waktu
subplot(3,1,2);
plot(time_dst, soc_dst, 'k');
ylabel('SOC (%)');
xlabel('Time (s)');
title('SOC');
grid on;

% Baris 3: Uocv terhadap waktu
subplot(3,1,3);
plot(time_dst, uoc_dst, 'm');
ylabel('U_{ocv} (V)');
xlabel('Time (s)');
title('Open Circuit Voltage (U_{ocv})');
grid on;

% Hitung error E(k) = Vt - Uoc
E = voltage_dst - uoc_dst;
N = length(E);

% Inisialisasi regresor Phi dan target E(k)
Phi = zeros(N, 5);   % [E(k-1), E(k-2), I(k), I(k-1), I(k-2)]
E_target = zeros(N, 1);  % Y(k) = E(k)

for k = 1:N
    if k >= 2
        E_k_1 = E(k-1); I_k_1 = current_dst(k-1);
    else
        E_k_1 = 0; I_k_1 = 0;
    end

    if k >= 3
        E_k_2 = E(k-2); I_k_2 = current_dst(k-2);
    else
        E_k_2 = 0; I_k_2 = 0;
    end

    I_k = current_dst(k);

    Phi(k, :) = [E_k_1, E_k_2, I_k, I_k_1, I_k_2];
    E_target(k) = E(k); 
end

% Inisialisasi awal
lambda_min = 0.98;
h = 0.7;
e_base = 0.05;

n_param = 5;
N = size(Phi,1);

theta = zeros(n_param, N);
theta(:,1) = zeros(n_param, 1);
P = 1000 * eye(n_param);
e_all = zeros(N,1);
lambda_all = zeros(N,1);
K_all = zeros(n_param, N);
y_hat_all = zeros(N, 1);

lambda = lambda_min + (1 - lambda_min) * h^0;

for k = 1:N
    phi_k = Phi(k,:)';
    E_k = E_target(k);
    y_hat = theta(:,k)' * phi_k;
    error_k = E_k - y_hat;

    K = P * phi_k / (lambda + phi_k' * P * phi_k);

    if k < N
        theta(:,k+1) = theta(:,k) + K * error_k;
    end

    P = (1/lambda) * (P - K * phi_k' * P);

    % Simpan variabel iteratif
    e_all(k) = error_k;
    lambda_all(k) = lambda;
    K_all(:,k) = K;
    y_hat_all(k) = y_hat;

    % Update lambda sesuai paper
    epsilon_k = round((error_k / e_base)^2);
    lambda = lambda_min + (1 - lambda_min) * h^epsilon_k;
end


%===========================================================================%
%========================Lavenberg-Marquadt=================================%
%=======================Plot ke-5Theta====================================%
N = size(theta, 2);

figure;
for i = 1:5
    subplot(5,1,i);                  % Buat 5 subplot vertikal
    plot(1:N, theta(i,:), 'LineWidth', 1.2);
    ylabel(['\theta_' num2str(i)]);
    if i == 1
        title('Evolusi Setiap Parameter \theta');
    end
    if i == 5
        xlabel('Sample / Time Index');
    end
    grid on;
end
T = 5;    % Sesuaikan sampling time Anda, misal T = 50
N = size(theta, 2);

% Alokasi array
R0_direct = zeros(1, N);
R1_direct = zeros(1, N);
R2_direct = zeros(1, N);
C1_direct = zeros(1, N);
C2_direct = zeros(1, N);

R0_LM = zeros(1, N);
R1_LM = zeros(1, N);
R2_LM = zeros(1, N);
C1_LM = zeros(1, N);
C2_LM = zeros(1, N);

% Parameter inisialisasi untuk LM
params0 = [0.01, 0.01, 0.01, 1000, 1000];
lb = [0, 0, 0, 0, 0];
ub = [inf, inf, inf, inf, inf];

for k = 1:N
    %=== Ambil theta tiap waktu ===%
    th = theta(:, k);

    %===== 1. Direct Calculation =====%
    theta1 = th(1); theta2 = th(2); theta3 = th(3); theta4 = th(4); theta5 = th(5);

    a = (theta4 - theta3 - theta5) / (1 + theta1 - theta2);
    b = T^2 * (1 + theta1 - theta2) / (4 * (1 - theta1 - theta2));
    c = T * (1 + theta2) / (1 - theta1 - theta2);
    d = -(theta3 - theta4 - theta5) / (1 - theta1 - theta2);
    f = T * (theta5 - theta3) / (1 - theta1 - theta2);

    tau1 = (c + sqrt(c^2 - 4*b)) / 2;
    tau2 = (c - sqrt(c^2 - 4*b)) / 2;
    R0_direct(k) = a;
    R1_direct(k) = (tau1 * (d - a) + a * c - f) / (tau1 - tau2);
    R2_direct(k) = d - a - R1_direct(k);
    C1_direct(k) = tau1 / R1_direct(k);
    C2_direct(k) = tau2 / R2_direct(k);

    %===== 2. LM (Levenberg-Marquardt) Fitting =====%
    try
        theta_exp_k = th;  % Target theta pada iterasi ke-k
        params_hat = lsqnonlin(@(p) rc_theta_error(p, theta_exp_k, T), params0, lb, ub, optimset('Display','off'));
        R0_LM(k) = params_hat(1);
        R1_LM(k) = params_hat(2);
        R2_LM(k) = params_hat(3);
        C1_LM(k) = params_hat(4);
        C2_LM(k) = params_hat(5);
    catch
        % Jika fitting gagal (misal, theta tidak valid), isi NaN
        R0_LM(k) = NaN; R1_LM(k) = NaN; R2_LM(k) = NaN; C1_LM(k) = NaN; C2_LM(k) = NaN;
    end
end

%--- Plot RC hasil Direct Calculation saja ---
figure('Name','Direct Calculation RC Tracking');
subplot(3,1,1);
plot(R0_direct, 'b', 'DisplayName','R0'); hold on;
plot(R1_direct, 'r', 'DisplayName','R1');
plot(R2_direct, 'g', 'DisplayName','R2');
legend; ylabel('\Omega'); title('Direct: Resistances'); grid on;

subplot(3,1,2);
plot(C1_direct, 'm', 'DisplayName','C1');
legend; ylabel('F'); title('Direct: C1'); grid on;

subplot(3,1,3);
plot(C2_direct, 'k', 'DisplayName','C2');
legend; ylabel('F'); title('Direct: C2'); xlabel('Sample'); grid on;


%--- Plot RC hasil LM saja ---
figure('Name','Levenberg-Marquardt RC Tracking (Per Parameter)');

subplot(5,1,1);
plot(R0_LM, 'b'); ylabel('R0 (\Omega)'); title('LM: R0'); grid on;

subplot(5,1,2);
plot(R1_LM, 'r'); ylabel('R1 (\Omega)'); title('LM: R1'); grid on;

subplot(5,1,3);
plot(R2_LM, 'g'); ylabel('R2 (\Omega)'); title('LM: R2'); grid on;

subplot(5,1,4);
plot(C1_LM, 'm'); ylabel('C1 (F)'); title('LM: C1'); grid on;

subplot(5,1,5);
plot(C2_LM, 'k'); ylabel('C2 (F)'); title('LM: C2'); grid on;
xlabel('Sample / Time Index');
start_idx = 1001;  % mulai plot dari iterasi ke-1001
end_idx = length(R0_LM);

figure('Name','Levenberg-Marquardt RC Tracking (R0, R1, R2, >1000 Iterasi)');

subplot(3,1,1);
plot(start_idx:end_idx, R0_LM(start_idx:end_idx), 'b');
ylabel('R0 (\Omega)'); title('LM: R0 (>1000 Iterasi)'); grid on;

subplot(3,1,2);
plot(start_idx:end_idx, R1_LM(start_idx:end_idx), 'r');
ylabel('R1 (\Omega)'); title('LM: R1 (>1000 Iterasi)'); grid on;

subplot(3,1,3);
plot(start_idx:end_idx, R2_LM(start_idx:end_idx), 'g');
ylabel('R2 (\Omega)'); title('LM: R2 (>1000 Iterasi)'); grid on;
xlabel('Sample / Time Index');


%======= Pastikan function rc2theta dan rc_theta_error sudah tersedia di file ======


%=============Persamaan (5) dari Paper================%
function theta_model = rc2theta(params, T)
    R0 = params(1); R1 = params(2); R2 = params(3);
    C1 = params(4); C2 = params(5);

    den = -T^2 - 2*T*(R1*C1 + R2*C2) - 4*R1*C1*R2*C2;
    theta1 = (2*T^2 - 8*R1*C1*R2*C2) / den;
    theta2 = (T^2 - 2*T*(R1*C1 + R2*C2) + 4*R1*C1*R2*C2) / den;
    theta3 = (T^2*(R0 + R1 + R2) + 2*T*(R0*R1*C1 + R0*R2*C2 + R2*R1*C2 + R2*R1*C1) + 4*R0*R1*C1*R2*C2) / den;
    theta4 = (2*T^2*(R0 + R1 + R2) - 8*R0*R1*C1*R2*C2) / den;
    theta5 = (T^2*(R0 + R1 + R2) - 2*T*(R0*R1*C1 + R0*R2*C2 + R1*R2*C2 + R2*R1*C1) + 4*R0*R1*C1*R2*C2) / den;

    theta_model = [theta1; theta2; theta3; theta4; theta5];
end

%============== Fungsi error untuk Lavenberg-Marquadt========%
function err = rc_theta_error(params, theta_exp, T)
    theta_model = rc2theta(params, T);
    err = theta_model - theta_exp;
end


function soc = calculate_soc(current, time, Q_nominal, soc_init)
    dt = [0; diff(time)];
    soc = zeros(length(current), 1);
    soc(1) = soc_init;

    for k = 2:length(current)
        delta_soc = (current(k) * dt(k)) / (Q_nominal * 3600) * 100;
        soc(k) = soc(k-1) + delta_soc;
    end

    soc = max(0, min(100, soc)); % Clamp SOC di antara 0–100%
end



