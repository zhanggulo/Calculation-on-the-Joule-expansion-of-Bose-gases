%% Bogoliubov-Popov Self-Consistent Calculation
%% Inversion Temperature
clear; clc

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

% 被积函数定义
int_func_1 = @(k, n0, a, t) ...
    k.^2 .* E_ki(k, n0, a) ./ (exp(E_ki(k, n0, a) ./ t) - 1);

int_func_2 = @(k, n0, a, t) ...
    k.^2 .* ((k.^2 / 2) ./ E_ki(k, n0, a)) .* ...
    (exp(E_ki(k, n0, a)./t) .* (1 - E_ki(k, n0, a)./t) - 1) ./ ...
    (exp(E_ki(k, n0, a) ./ t) - 1).^2;

% 主方程
equ = @(vars, a) [
    % 方程 1: 反转温度方程
    2 * pi * a + 2 * pi * a * nex(vars(1), a, vars(3))^2 ...
    + (128*sqrt(pi)/5) * a^(5/2) * vars(1)^(5/2) ...
    - quadgk(@(k) int_func_1(k, vars(1), a, vars(3)), k_min, k_max, 'RelTol', 1e-8) / (2 * pi^2) ...
    + (2 * a * vars(1) / pi) * quadgk(@(k) int_func_2(k, vars(1), a, vars(3)), k_min, k_max, 'RelTol', 1e-8);

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

%% Isoenergy Curves (fixed Ti)
clear; clc;

% 定义常量
T_c = 2 * pi / zeta(3/2)^(2/3); % 临界温度
Ti_values = (0.1:0.1:0.4) * T_c; % 不同的初始 Ti 值
x_values = linspace(1, 8, 100); % x = V_f / V_i 的范围
Tf = zeros(length(Ti_values), length(x_values)); % 存储结果
Tf_over_Tc_results = zeros(length(Ti_values), length(x_values)); % 存储结果
n0_values = zeros(length(Ti_values), length(x_values)); % 存储每条线的 n0
nex_values = zeros(length(Ti_values), length(x_values)); % 存储每条线的 nex
k_min = 0; % 积分下限
k_max = 10; % 积分上限
a = 0.03; % 固定的 a 值

% 色散关系
E_ki = @(k, n0, a) sqrt((k.^2 / 2) .* (k.^2 / 2 + 8 * pi * a * n0)); % 初态色散关系
E_kf = @(k, n0, a, x) sqrt((k.^2 / 2) .* (k.^2 / 2 + 8 * pi * a * n0 / x)); % 终态色散关系

% 定义 nex
nex = @(n0, a, t) quadgk(@(k) k.^2 .* ...
    ((k.^2 / 2 + 4 * pi * a * n0) ./ E_ki(k, n0, a)) ./ ...
    (exp(E_ki(k, n0, a) ./ t) - 1), k_min, k_max, 'RelTol', 1e-8) / (2 * pi^2);

% 被积函数定义
int_func_ki = @(k, n0, a, Ti) ...
    k.^2 .* E_ki(k, n0, a) ./ (exp(E_ki(k, n0, a) ./ Ti) - 1);
int_func_kf = @(k, n0, a, Tf, x) ...
    k.^2 .* E_kf(k, n0, a, x) ./ (exp(E_kf(k, n0, a, x) ./ Tf) - 1);

% 主方程
equ = @(vars, x, Ti) [
    % 方程 1: 等能条件
    quadgk(@(k) int_func_kf(k, vars(1), a, vars(3), x), k_min, k_max, 'RelTol', 1e-8) ...
    - (1 / x) * quadgk(@(k) int_func_ki(k, vars(1), a, Ti), k_min, k_max, 'RelTol', 1e-8) ...
    - 4 * pi^3 * a * ((1 / x) * (1 - 1 / x) + nex(vars(1), a, vars(3))^2 * (1 / x) * (1 - 1 / x) ...
    + 128 / (15 * sqrt(pi)) * a^(3/2) * vars(1)^(5/2) * (1 / x) * (1 - (1 / x)^(3/2)));

    % 方程 2: 自洽条件 (n0 + nex = 1)
    vars(1) + nex(vars(1), a, vars(3)) - 1;

    % 方程 3: nex 的定义
    vars(2) - nex(vars(1), a, vars(3));
];

% 初始猜测值
vars_initial = [0.5, 0.5, 0.5 * T_c]; % [n0_guess, nex_guess, Tf_guess]

% 数值求解
for j = 1:length(Ti_values)
    Ti = Ti_values(j); % 当前 Ti
    for i = 1:length(x_values)
        x = x_values(i); % 当前 x
        try
            % 用 fsolve 求解方程
            vars_solution = fsolve(@(vars) equ(vars, x, Ti), vars_initial, ...
                                    optimoptions('fsolve', 'Display', 'off'));
            n0_values(j, i) = vars_solution(1);
            nex_values(j, i) = vars_solution(2);
            Tf(j, i) = vars_solution(3);
            Tf_over_Tc_results(j, i) = Tf(j, i) / T_c;
            Tf_over_T_results(j, i) = Tf(j, i) / Ti;
            fprintf('Ti/Tc = %.1f, x = %.2f, n0 = %.4f, nex = %.4f, Tf/Tc = %.4f\n', ...
                    Ti / T_c, x, vars_solution(1), vars_solution(2), vars_solution(3) / T_c);
        catch ME
            warning('Ti/Tc = %.1f, x = %.2f, 求解失败: %s', Ti / T_c, x, ME.message);
            n0_values(j, i) = NaN;
            nex_values(j, i) = NaN;
            Tf_over_Tc_results(j, i) = NaN;
        end
    end
end

% 绘图
figure;
hold on;
colors = lines(length(Ti_values));
for j = 1:length(Ti_values)
    plot(x_values, Tf_over_T_results(j, :), 'LineWidth', 1, 'Color', colors(j, :), ...
         'DisplayName', sprintf('$T_i / T_c = %.1f$', Ti_values(j) / T_c));
end
xlabel('$V_f / V_i$', 'Interpreter', 'latex');
ylabel('$T_f / T_{c,i}$', 'Interpreter', 'latex');
hold off;
box on;

%% %% Isoenergy Curves (fixed Vf / Vi)
clear; clc;

% 定义常量
T_c = 2 * pi / zeta(3/2)^(2/3); % 临界温度
Ti_values = linspace(0.0, 0.5, 100) * T_c; % Ti 的范围
x = 2.6; % 固定膨胀比 V_f / V_i
Tf_results = zeros(size(Ti_values)); % 存储 Tf 结果
Tf_over_Ti_results = zeros(size(Ti_values)); % 存储 Tf / Ti 结果
n0_values = zeros(size(Ti_values)); % 存储每个 Ti 的 n0
nex_values = zeros(size(Ti_values)); % 存储每个 Ti 的 nex
k_min = 0; % 积分下限
k_max = 10; % 积分上限
a = 0.04; % 固定的 a 值

% 色散关系
E_ki = @(k, n0, a) sqrt((k.^2 / 2) .* (k.^2 / 2 + 8 * pi * a * n0)); % 初态色散关系
E_kf = @(k, n0, a, x) sqrt((k.^2 / 2) .* (k.^2 / 2 + 8 * pi * a * n0 / x)); % 终态色散关系

% 定义 nex
nex = @(n0, a, t) quadgk(@(k) k.^2 .* ...
    ((k.^2 / 2 + 4 * pi * a * n0) ./ E_ki(k, n0, a)) ./ ...
    (exp(E_ki(k, n0, a) ./ t) - 1), k_min, k_max, 'RelTol', 1e-8) / (2 * pi^2);

% 被积函数定义
int_func_ki = @(k, n0, a, Ti) ...
    k.^2 .* E_ki(k, n0, a) ./ (exp(E_ki(k, n0, a) ./ Ti) - 1);
int_func_kf = @(k, n0, a, Tf, x) ...
    k.^2 .* E_kf(k, n0, a, x) ./ (exp(E_kf(k, n0, a, x) ./ Tf) - 1);

% 主方程
equ = @(vars, x, Ti) [
    % 方程 1: 等能条件
    quadgk(@(k) int_func_kf(k, vars(1), a, vars(3), x), k_min, k_max, 'RelTol', 1e-8) ...
    - (1 / x) * quadgk(@(k) int_func_ki(k, vars(1), a, Ti), k_min, k_max, 'RelTol', 1e-8) ...
    - 4 * pi^3 * a * ((1 / x) * (1 - 1 / x) + nex(vars(1), a, vars(3))^2 * (1 / x) * (1 - 1 / x) ...
    + 128 / (15 * sqrt(pi)) * a^(3/2) * vars(1)^(5/2) * (1 / x) * (1 - (1 / x)^(3/2)));

    % 方程 2: 自洽条件 (n0 + nex = 1)
    vars(1) + nex(vars(1), a, vars(3)) - 1;

    % 方程 3: nex 的定义
    vars(2) - nex(vars(1), a, vars(3));
];

% 初始猜测值
vars_initial = [0.5, 0.5, 0.5 * T_c]; % [n0_guess, nex_guess, Tf_guess]

% 数值求解
for i = 1:length(Ti_values)
    Ti = Ti_values(i); % 当前 Ti
    try
        % 用 fsolve 求解方程
        vars_solution = fsolve(@(vars) equ(vars, x, Ti), vars_initial, ...
                                optimoptions('fsolve', 'Display', 'off'));
        n0_values(i) = vars_solution(1);
        nex_values(i) = vars_solution(2);
        Tf_results(i) = vars_solution(3);
        Tf_over_Ti_results(i) = vars_solution(3) / Ti;
        fprintf('Ti/Tc = %.4f, n0 = %.4f, nex = %.4f, Tf/Ti = %.4f\n', ...
                Ti / T_c, vars_solution(1), vars_solution(2), vars_solution(3) / Ti);
    catch ME
        warning('Ti/Tc = %.4f, 求解失败: %s', Ti / T_c, ME.message);
        n0_values(i) = NaN;
        nex_values(i) = NaN;
        Tf_over_Ti_results(i) = NaN;
    end
end

% 绘图
figure;
plot(Ti_values / T_c, Tf_over_Ti_results, 'LineWidth', 1);
xlabel('$T_i / T_{c,i}$', 'Interpreter', 'latex');
ylabel('$T_f / T_i$', 'Interpreter', 'latex');
grid on;
box on;

figure;
% 画出原始曲线
plot(Ti_values / T_c, Tf_results / T_c, 'LineWidth', 1); hold on;
plot(Ti_values / T_c, Ti_values / T_c, '--k', 'LineWidth', 1);

% 填充红色区域（虚线到实线右侧区域）
x_fill_red = Ti_values(Tf_results > Ti_values) / T_c; % x 轴值
y_fill_red = Tf_results(Tf_results > Ti_values) / T_c; % 实线值
fill([x_fill_red, fliplr(x_fill_red)], ...
     [x_fill_red, fliplr(y_fill_red)], 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% 填充蓝色区域（虚线到实线左侧区域）
x_fill_blue = Ti_values(Tf_results < Ti_values) / T_c; % x 轴值
y_fill_blue = Tf_results(Tf_results < Ti_values) / T_c; % 实线值
fill([x_fill_blue, fliplr(x_fill_blue)], ...
     [y_fill_blue, fliplr(x_fill_blue)], 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

hold off;
xlabel('$T_i / T_{c,i}$', 'Interpreter', 'latex');
ylabel('$T_f / T_{c,i}$', 'Interpreter', 'latex');
box on;