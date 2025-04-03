clear; clc; close all;

a = 0.01; % 密度参数
k_max = 10; % 截断积分上限，避免发散
T_c = 2 * pi / zeta(3/2)^(2/3); % 临界温度
t = linspace(0, 1, 100) * T_c; % 温度范围 [0, T_c]

E_k = @(k, a) sqrt((k.^2 / 2) .* (k.^2 / 2 + 8 * pi * a));

% 被积函数
int_func = @(k, a, t) k.^2 .* E_k(k, a) ./ (exp(E_k(k, a) ./ t) - 1);

% 数值积分项定义（避免无穷积分）
int_term = @(t, a) arrayfun(@(temp) ...
    integral(@(k) int_func(k, a, temp), 0, k_max, ...
    'AbsTol', 1e-6, 'RelTol', 1e-4), t);

% 总能量数值计算
E = @(t, a) 2 .* pi .* a .* (1 + 128./(15.*sqrt(pi)) .* a.^(3/2)) + ...
    int_term(t, a) ./ (2 * pi)^2;

% 近似解析解（低温极限）
E1 = @(t, a) 2 .* pi .* a .* (1 + 128./(15.*sqrt(pi)) .* a.^(3/2)) + ...
    pi^2 .* t.^4 ./ (60*(4*pi*a)^(3/2));

E_vals = E(t, a); % 数值解
E1_vals = E1(t, a); % 近似解

figure;
plot(t / T_c, E_vals, 'LineWidth', 1); hold on;
plot(t / T_c, E1_vals, 'LineWidth', 1);
legend('Numerical', 'Analytical', 'Interpreter', 'latex');
xlabel('$T / T_c$', 'Interpreter', 'latex');
ylabel('$E$', 'Interpreter', 'latex');
title('$E(n=0.01,T)$', 'Interpreter', 'latex');
grid on;
box on;

%% vary n
clear; clc; close all;

a_values = linspace(0, 0.05, 100); % 密度范围 [0, 0.05]
k_max = inf; % 截断积分上限，避免发散
T_c = 2 * pi / zeta(3/2)^(2/3); % 临界温度
t = 0.1; % 固定温度 (单位：T_c)

E_k = @(k, a) sqrt((k.^2 / 2) .* (k.^2 / 2 + 8 * pi * a));

% 分布函数
f = @(k, a, t) 1 ./ (exp(E_k(k, a) ./ t) - 1);

% 被积函数
int_func = @(k, a, t) k.^2 .* E_k(k, a) .* f(k, a, t);

% 数值积分项定义（避免无穷积分）
int_term = @(t, a) integral(@(k) int_func(k, a, t), 0, k_max, ...
    'AbsTol', 1e-6, 'RelTol', 1e-4);

% 总能量数值计算
E = @(t, a) 2 .* pi .* a .* (1 + 128./(15.*sqrt(pi)) .* a.^(3/2)) + ...
    int_term(t, a) ./ (2 * pi)^2;

% 近似解析解（低温极限）
E1 = @(t, a) 2 .* pi .* a .* (1 + 128./(15.*sqrt(pi)) .* a.^(3/2)) + ...
    pi^2 .* t.^4 ./ (60*(4*pi*a)^(3/2));

% 近似解析解（高温极限）
E2 = @(t, a) 2 .* pi .* a .* (1 + 128./(15.*sqrt(pi)) .* a.^(3/2)) + ...
    t * gamma(5/2) * zeta(5/2) / (2 * pi^2);

E_vals = arrayfun(@(a) E(t * T_c, a), a_values); % 数值解
E1_vals = arrayfun(@(a) E1(t * T_c, a), a_values); % 近似解
E2_vals = arrayfun(@(a) E2(t * T_c, a), a_values);

figure;
plot(a_values, E_vals, 'LineWidth', 1); hold on;
plot(a_values, E1_vals, '--', 'LineWidth', 1);
plot(a_values, E2_vals, '-.', 'LineWidth', 1);

legend('Numerical', 'Analytical (small k)', 'Analytical (large k)');
xlabel('$n^{1/3}a$', 'Interpreter', 'latex');
ylabel('$E(a, T)/(\hbar^2 n^{2/3}/m)$', 'Interpreter', 'latex');
title('$E(a, T = 0.1T_c)$', 'Interpreter', 'latex');
grid on;
box on;