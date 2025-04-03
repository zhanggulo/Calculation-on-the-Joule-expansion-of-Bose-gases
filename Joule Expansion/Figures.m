%% Inversion Temperature
clear; clc

eta = 1; % Condensation Fraction
x = linspace(0, 0.05, 1000); % Initial density in unit of (1/a)^3

% Perturbation
tI_P = ((2 - eta^2) * 2 * zeta(3/2) ^ (5/3) / (3 * zeta(5/2)) * x) .^ (2/5);

% Bogoliubov (Small Momentum Approximation)
tI_BS = ((3/4)*((4*pi)^(5/2)/pi^6)*zeta(3/2)^(8/3) .* x.^(5/2) .* (1/2 + 32/(5*sqrt(pi)) .* x.^(3/2))) .^ (1/4);

% Bogoliubov (Large Momentum Approximation)
tI_BL = (sqrt(pi/2) * (2*zeta(3/2)^(5/3)/(gamma(5/2)*zeta(5/2))) * x .* (1 + 64/(5*sqrt(pi)) * x.^(3/2))) .^ (2/5);

figure;
hold on;
plot(x, tI_BS, "LineWidth", 1);
plot(x, 2*zeta(3/2)^(2/3)*x, 'LineWidth', 1);
hold off;
xlabel('$n^{1/3}a$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('$T_I/T_c$', 'Interpreter', 'latex')

%% Isoenergy Curves (Ti/Tc was fixed)
clear; clc

eta = 1; % Condensation Fraction
x = linspace(1, 8, 100); % Expansion fraction Vf / Vi
ni = 0.03; % Initial density in unit of (1/a)^3
ti = 0.1; % Initial temperature in unit of Tc

% Perturbation
tf = @(ti) ((ti) .^ (5/2) ./ x + 2 .* zeta(3/2) .^ (5/3) / (3 * zeta(5/2)) .* ...
    (2 - eta.^2) .* ni ./ x .* (1 - 1 ./ x)) .^ (2/5); % Equation of equal energy

% figure % Perturbation
hold on;
for i = 1 : length(ti)
plot(x, tf(ti(i))/ti(i),'r--', "LineWidth", 0.6);
end
hold off
xlabel('$V_f/V_i$', 'Interpreter','latex')
ylabel('$T_f/T_{c,i}$','Interpreter','latex')

%% Isoenergy Curves (Vi/Vf was fixed)
clear; clc

% 参数定义
eta = 1; % 冷凝分数 (Condensation Fraction)
x_vals = 5; % 膨胀比 Vf / Vi
ni = 0.03; % 初始密度，单位为 (1/a)^3
ti = linspace(0.1, 0.5, 1000); % 初始温度 Ti，单位为 Tc

% 计算 tf_P / ti
% figure; 
hold on;
colors = lines(length(x_vals)); % 生成颜色
for i = 1:length(x_vals)
    x = x_vals(i);
    tf = ((ti) .^ (5/2) ./ x + 2 * zeta(3/2) .^ (5/3) / (3 * zeta(5/2)) .* ...
        (2 - eta^2) .* ni ./ x .* (1 - 1 ./ x)) .^ (2/5); % 计算 Tf

    % 绘图
    plot(ti, tf ./ ti,'r--', 'LineWidth', 0.6);
end

% 轴标签和标题
xlabel('$T_i/T_{c,i}$', 'Interpreter', 'latex', 'FontSize', 10);
ylabel('$T_f/T_i$', 'Interpreter', 'latex', 'FontSize', 10);

% 格式美化
box on;
set(gca, 'FontSize', 10);
hold off;