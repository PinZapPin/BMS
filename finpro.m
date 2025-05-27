clc;
clearvars;

% Mencari OCV vs SOC
data = readtable('007-OCV.csv');
Q_nominal = 1.1;
time = data.Step_Time_s;
current = data.Current_A;
voltage = data.Voltage_V;
dt = [0; diff(time)];

idx_cha = current > 0;
idx_dis = current < 0;

current_cha = current(idx_cha);
voltage_cha = voltage(idx_cha);
time_cha = time(idx_cha);
dt_cha = dt(idx_cha);

current_dis = current(idx_dis);
voltage_dis = voltage(idx_dis);
time_dis = time(idx_dis);
dt_dis = dt(idx_dis);

soc_cha = zeros(length(current_cha),1);
soc_dis = zeros(length(current_dis),1);

for i = 2:length(current_cha)
    soc_cha(i) = soc_cha(i-1) + (current_cha(i) * dt_cha(i) / (Q_nominal * 3600)) * 100;
end

soc_dis(1) = 100;
for i = 2:length(current_dis)
    soc_dis(i) = soc_dis(i-1) + (current_dis(i) * dt_dis(i) / (Q_nominal * 3600)) * 100;
end

soc_dis_norm = soc_dis / 100;
p_dis = polyfit(soc_dis_norm, voltage_dis, 9);
soc_dis_fit = linspace(min(soc_dis_norm), max(soc_dis_norm), 200);
voltage_dis_fit = polyval(p_dis, soc_dis_fit);

soc_cha_norm = soc_cha / 100;
p_cha = polyfit(soc_cha_norm, voltage_cha, 9);
soc_cha_fit = linspace(min(soc_cha_norm), max(soc_cha_norm), 200);
voltage_cha_fit = polyval(p_cha, soc_cha_fit);

subplot(1,2,1);
scatter(soc_cha_norm*100, voltage_cha, 10, 'r', 'filled');
hold on;
plot(soc_cha_fit*100, voltage_cha_fit, 'b-', 'LineWidth', 2);
xlabel('SOC (%)');
ylabel('Voltage (V)');
title('Charge: Polynomial Fit vs Measured Voltage');
grid on;
hold off;

subplot(1,2,2);
scatter(soc_dis_norm*100, voltage_dis, 10, 'r', 'filled');
hold on;
plot(soc_dis_fit*100, voltage_dis_fit, 'b-', 'LineWidth', 2);
xlabel('SOC (%)');
ylabel('Voltage (V)');
title('Discharge: Polynomial Fit vs Measured Voltage');
grid on;
hold off;

