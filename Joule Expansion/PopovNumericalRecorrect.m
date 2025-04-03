%% Inversion Temperature (Modified)
clear; clc;

% 定义常量
T_c = 2 * pi / zeta(3/2)^(2/3); % 临界温度
a_values = linspace(0, 0.05, 1000); % 密度
t_values = zeros(size(a_values)); % 存储每个 a 对应的温度
n0_values = zeros(size(a_values)); % 存储每个 a 对应的 n0
nex_values = zeros(size(a_values)); % 存储每个 a 对应的 nex
k_min = 0; % 积分下限
k_max = 10; % 积分上限

% 色散关系
E_ki = @(k, n0, a) sqrt((k.^2 / 2) .* (k.^2 / 2 + 8 * pi * a * n0));

% 定义 nex
nex = @(n0, a, t) quadgk(@(k) k.^2 .* ...
    ((k.^2 / 2 + 4 * pi * a * n0) ./ E_ki(k, n0, a)) ./ ...
    (exp(E_ki(k, n0, a) ./ t) - 1), k_min, k_max, 'RelTol', 1e-8) / (2 * pi^2);

% 计算 \partial n_0 / \partial V
partial_n0_V = @(n0, a, t) -1 ./ ( (1 / (2 * pi^2)) * quadgk(@(k) k.^2 .* ...
    ((4 * pi * a .* E_ki(k, n0, a) - (k.^2 / 2 + 4 * pi * a * n0) .* (4 * pi * a .* k.^2 / 2) ./ E_ki(k, n0, a)) ./ ...
    (E_ki(k, n0, a).^2 .* (exp(E_ki(k, n0, a) ./ t) - 1)) - ...
    ((k.^2 / 2 + 4 * pi * a * n0) .* exp(E_ki(k, n0, a) ./ t) .* (4 * pi * a .* k.^2 / 2) ./ (t * E_ki(k, n0, a).^2)) ./ ...
    (exp(E_ki(k, n0, a) ./ t) - 1).^2), k_min, k_max, 'RelTol', 1e-8) + 1);

% 计算 \partial n_{ex} / \partial V
partial_nex_V = @(n0_V) -1 - n0_V;

% 被积函数定义
int_func_1 = @(k, n0, a, t) ...
    k.^2 .* E_ki(k, n0, a) ./ (exp(E_ki(k, n0, a) ./ t) - 1);

int_func_2 = @(k, n0, a, t, dn0_dV) ...
    k.^2 .* dn0_dV .* (4 * pi * a .* k.^2 / 2) ./ E_ki(k, n0, a) .* ...
    (exp(E_ki(k, n0, a)./t) .* (1 - E_ki(k, n0, a)./t) - 1) ./ ...
    (exp(E_ki(k, n0, a) ./ t) - 1).^2;

% 方程 1: 主方程
equ = @(vars, a) [
    -2 * pi * a + 2 * pi * a * nex(vars(1), a, vars(3))^2 ...
    + 4 * pi * a * nex(vars(1), a, vars(3)) * partial_nex_V(partial_n0_V(vars(1), a, vars(3))) ...
    + 2 * pi * a * (128 / sqrt(pi)) * a^(3/2) * vars(1)^(5/2) ...
    + 4 * pi * a^(5/2) * (32 / 3) * vars(1)^(3/2) * partial_n0_V(vars(1), a, vars(3)) ...
    + quadgk(@(k) int_func_1(k, vars(1), a, vars(3)), k_min, k_max, 'RelTol', 1e-8) / (2 * pi^2) ...
    + quadgk(@(k) int_func_2(k, vars(1), a, vars(3), partial_n0_V(vars(1), a, vars(3))), k_min, k_max, 'RelTol', 1e-8) / (2 * pi^2);

    % 方程 2: 自洽条件 (n0 + nex = 1)
    vars(1) + nex(vars(1), a, vars(3)) - 1;

    % 方程 3: nex 的定义
    vars(2) - nex(vars(1), a, vars(3));
    ];

% 初始猜测值
vars_initial = [0.5, 0.5, 0.5 * T_c]; % [n0_guess, nex_guess, t_guess]

% 数值求解
for i = 1:length(a_values)
    a = a_values(i); % 当前密度
    if a == 0 % 单独设定 a=0 的情形 (可以解析求解)
        n0_values(i) = 1;
        nex_values(i) = 0;
        t_values(i) = 0;
        continue;
    end
    try
        vars_solution = fsolve(@(vars) equ(vars, a), vars_initial, optimoptions('fsolve', 'Display', 'off'));
        n0_values(i) = vars_solution(1);
        nex_values(i) = vars_solution(2);
        t_values(i) = vars_solution(3);
        fprintf('a = %.6f, n0 = %.4f, nex = %.4f, t = %.4f\n', a, vars_solution(1), vars_solution(2), vars_solution(3)/T_c);
    catch ME
        warning('a = %.6f, 求解失败: %s', a, ME.message);
        n0_values(i) = NaN;
        nex_values(i) = NaN;
        t_values(i) = NaN;
    end
end

T_over_Tc = t_values / T_c; % 归一化温度单位

% 绘图
figure;
plot(a_values, T_over_Tc, 'LineWidth', 1);
xlabel('$n^{1/3}a$', 'Interpreter', 'latex');
ylabel('$T_{I} / T_c$', 'Interpreter', 'latex');
title('Inversion Temperature');
box on;

figure;
yyaxis left;
plot(a_values, n0_values, 'LineWidth', 1);
ylabel('$n_0/n$', 'Interpreter', 'latex');
yyaxis right
plot(a_values, nex_values, 'LineWidth', 1);
ylabel('$\tilde{n}/n$', 'Interpreter', 'latex');
xlabel('$n^{1/3}a$', 'Interpreter', 'latex');
box on;
