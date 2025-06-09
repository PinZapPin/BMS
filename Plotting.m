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

% Plot 1: Tegangan dan Arus terhadap waktu
figure;
yyaxis left;
plot(time_dst, voltage_dst, 'b');
ylabel('Voltage (V)');
yyaxis right;
plot(time_dst, current_dst, 'r');
ylabel('Current (A)');
xlabel('Time (s)');
title('Voltage and Current');
grid on;

% Plot 2: SOC terhadap waktu
figure;
plot(time_dst, soc_dst, 'k');
ylabel('SOC (%)');
xlabel('Time (s)');
title('State of Charge (SOC)');
grid on;

% Plot 3: Open Circuit Voltage terhadap waktu
figure;
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
h = 0.9;
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

T = 5;

% --- Parameter LM ---
params_init = [0, 0, 0, 0, 0]; % Cek kisaran yang kira-kira mungkin
lb = [1e-7, 1e-7, 1e-7, 10, 10];
ub = [10,10, 10, 2e5, 2e5];
mu_init = 1e-4;       % damping LM
max_iter = 50;        % iterasi maksimum LM per titik
tol = 1e-6;           % toleransi konvergensi
delta = 1e-3;         % langkah Jacobian numerik

% --- Untuk hasil tracking ---
RC_tracking = nan(N, 5); % hasil LM tiap sample
err_tracking = nan(N, 1); % error fitting LM tiap sample

%==========Loop Theta dan Parameter RC menggunakan LM ===========%
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

    % --- LM Start (untuk waktu k) ---
    theta_exp = theta(:,k);   
    params = params_init;
    mu = mu_init;

    for iter = 1:max_iter
        theta_model = rc2theta(params, T);
        err = theta_model - theta_exp;

        % Hitung Jacobian
        J = zeros(5,5);
        for j = 1:5
            params_temp = params;
            params_temp(j) = params_temp(j) + delta;
            theta_temp = rc2theta(params_temp, T);
            J(:,j) = (theta_temp - theta_model) / delta;
        end

        % LM update
        H = J'*J + mu*eye(5);         % Hessian semu
        dp = H \ (J'*err);            % langkah update (delta param)

        params_new = params - dp';     % update param
        % Boundary check:
        params_new = max(params_new, lb);
        params_new = min(params_new, ub);

        theta_new = rc2theta(params_new, T);
        err_new = theta_new - theta_exp;

        if norm(err_new) < norm(err)
            params = params_new;
            mu = mu / 10;
            if norm(err - err_new) < tol
                break;  % stop iterasi jika sudah cukup konvergen
            end
        else
            mu = mu * 10;
        end

    end
    
    % ---- Di sini tambahkan penyimpanan hasil params ke RC_tracking ----
    RC_tracking(k,:) = params;
    err_tracking(k) = norm(rc2theta(params, T) - theta_exp);

end

% ===================== Plot Lambda Tracking ===================== %
figure;
set(gcf, 'Position', [100 100 1000 400]);  % opsional: perbesar tampilan
plot(1:N, lambda_all, 'LineWidth', 1.5);
xlabel('Sample / Time Index');
ylabel('\lambda_k');
title('Evolusi Forgetting Factor \lambda_k');
grid on;


%=======================Fitting Theta terakhir via LM===============
N = size(theta,2);             % Pastikan N sesuai dengan theta hasil RLS/FFRLS
theta_exp2 = theta(:,N);       % Gunakan theta terakhir

err_hist2 = zeros(max_iter,1);
params2 = params_init;
mu2 = mu_init;

for iter2 = 1:max_iter
    theta_model2 = rc2theta(params2, T);
    err2 = theta_model2 - theta_exp2;
    err_hist2(iter2) = norm(err2);

    % Hitung Jacobian numerik
    J2 = zeros(5,5);
    for j2 = 1:5
        params_temp2 = params2;
        params_temp2(j2) = params_temp2(j2) + delta;
        theta_temp2 = rc2theta(params_temp2, T);
        J2(:,j2) = (theta_temp2 - theta_model2) / delta;
    end

    % LM update
    H2 = J2'*J2 + mu2*eye(5);
    dp2 = H2 \ (J2'*err2);

    params_new2 = params2 - dp2';   % transpose agar baris
    % Boundary check
    params_new2 = max(params_new2, lb);
    params_new2 = min(params_new2, ub);

    theta_new2 = rc2theta(params_new2, T);
    err_new2 = theta_new2 - theta_exp2;

    if norm(err_new2) < norm(err2)
        params2 = params_new2;
        mu2 = mu2 / 10;
        if norm(err2 - err_new2) < tol
            err_hist2 = err_hist2(1:iter2); % potong sisa zero jika konvergen duluan
            break;
        end
    else
        mu2 = mu2 * 10;
    end
end

RC_tracking2 = params2;
err_tracking2 = norm(rc2theta(params2, T) - theta_exp2);

% --- Display hasil parameter RC
disp('Parameter RC hasil LM (theta terakhir):');
fprintf('R0 = %.6f Ohm\n', params2(1));
fprintf('R1 = %.6f Ohm\n', params2(2));
fprintf('R2 = %.6f Ohm\n', params2(3));
fprintf('C1 = %.6f F\n',   params2(4));
fprintf('C2 = %.6f F\n',   params2(5));

% --- Plot error konvergensi LM
figure;
plot(err_hist2, 'LineWidth', 1.5);
xlabel('Iterasi LM');
ylabel('Norm Error');
title('Konvergensi Error LM (Fitting theta terakhir)');
grid on;

%================== Plot RC_tracking per parameter (1 grafik per jendela) ===============%
param_names = {'R0', 'R1', 'R2', 'C1', 'C2'};

for i = 1:5
    figure;
    set(gcf, 'Position', [100 100 1000 400]);  % opsional: perbesar tampilan
    plot(RC_tracking(:,i), 'LineWidth', 1.5);
    ylabel(param_names{i});
    xlabel('Sample / Time Index');
    title(['Tracking Parameter ' param_names{i}]);
    grid on;
end

figure;
plot(err_tracking, 'LineWidth', 1.5); grid on;
xlabel('Sample'); ylabel('Norm Fitting Error');
title('LM Error per Sample');



%=======================Plot ke-5 Theta secara terpisah ====================%
N = size(theta, 2);

for i = 1:5
    figure;
    set(gcf, 'Position', [100 100 1000 400]);  % (opsional) perbesar ukuran figure
    plot(1:N, theta(i,:), 'LineWidth', 1.5);
    ylabel(['\theta_' num2str(i)]);
    xlabel('Sample / Time Index');
    title(['Evolusi Parameter \theta_' num2str(i)]);
    grid on;
end

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



