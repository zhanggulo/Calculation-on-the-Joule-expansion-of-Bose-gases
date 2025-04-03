%% Invertion Temperature (Numerical Bogoliubov)
clear; clc;

% 常数定义
a_values = linspace(0, 0.05, 1000); % 初始密度
t_values = zeros(size(a_values)); % 存储每个 a 对应的 t
t_over_tc_values = zeros(size(a_values)); % 存储每个 a 对应的 T/Tc
T_c = 2 * pi / zeta(3/2)^(2/3); % 临界温度
k_min = 0; % 积分下限
k_max = 10; % 积分上限（避免无穷发散）

% 色散关系
E_k = @(k, a) sqrt((k.^2 / 2) .* (k.^2 / 2 + 8 * pi * a));

% 被积函数1
int_func_1 = @(k, a, t) k.^2 .* E_k(k, a) ./ (exp(E_k(k, a) ./ t) - 1);

% 被积函数2
int_func_2 = @(k, a, t) k.^2 .* ((k.^2 / 2) ./ E_k(k, a)) .* ...
    (exp(E_k(k, a)./t) .* (1 - E_k(k, a)./t) - 1) ./ (exp(E_k(k, a) ./ t) - 1).^2;

% 数值积分
int_term_1 = @(t, a) quadgk(@(k) int_func_1(k, a, t), k_min, k_max, 'RelTol', 1e-8);
int_term_2 = @(t, a) quadgk(@(k) int_func_2(k, a, t), k_min, k_max, 'RelTol', 1e-8);

% 定义反转温度方程方程
equ = @(t, a) 2 * pi * a + (128 * sqrt(pi) / 5) * a^(5/2) - ...
    int_term_1(t, a) / (2*pi^2) + (2*a/pi) * int_term_2(t, a);

% 求解温度 t 对每个 a
t_initial = 0.5 * T_c; % 初始猜测值
options = optimset('Display', 'off', 'TolX', 1e-6);

for i = 1:length(a_values)
    a = a_values(i); % 当前密度
    if a == 0 % 单独设定 a=0 的情形 (可以解析求解)
        t_values(i) = 0;
        continue;
    end
    try
        t_solution = fzero(@(t) equ(t, a), t_initial, options);
        t_values(i) = t_solution; % 存储温度
        t_over_tc_values(i) = t_solution / T_c; % 转换为 T / Tc
        fprintf('a = %.6f, t = %.4f, T/Tc = %.4f\n', a, t_solution, t_solution / T_c);
    catch ME
        warning('a = %.6f, 求解失败: %s', a, ME.message);
        t_values(i) = NaN; % 如果求解失败，记录为 NaN
        t_over_tc_values(i) = NaN;
    end
end

y = @(x) 2 * zeta(3/2)^(2/3) * x; % 声子谱近似条件

% 绘制结果
figure;
plot(a_values, t_over_tc_values, 'LineWidth', 1); hold on;
plot(a_values, y(a_values), '--k', 'LineWidth', 1);
xlabel('$n^{1/3}a$', 'Interpreter', 'latex');
ylabel('$T_{I} / T_c$', 'Interpreter', 'latex');
title('Bogoliubov Inversion Temperature (Numerical)');
box on;

%% Isoenergy Curves (Numerical Bogoliubov)
clear; clc;

% 基础参数
x = linspace(1, 5, 100); % Expansion fraction Vf / Vi
a = 0.05; % Initial density in unit of 1/n_i^(1/3)
T_c_i = 2 * pi / zeta(3/2)^(2/3); % 初始临界温度
ti_values = (0.35 : 0.05 : 0.5) * T_c_i; % 初始温度的不同取值
k_max = Inf; % 数值积分的上限 

% 初始化存储变量
tf_results = zeros(length(ti_values), length(x));

% 外部循环：不同 ti 值
for j = 1:length(ti_values)
    ti = ti_values(j);
    
    % 初始温度积分
    E_k_i = @(k) sqrt(k.^2/2 .* (k.^2/2 + 8*pi*a));
    int_ti_func = @(k) k.^2 .* E_k_i(k) ./ (exp(E_k_i(k) ./ ti) - 1);
    int_ti_func_term = integral(int_ti_func, 0, k_max, 'ArrayValued', true);
    
    % 初始化每个 ti 下的变量
    for i = 1:length(x)
        % 存储每个 x(i) 的终函数
        E_k_f = @(k) sqrt(k.^2/2 .* (k.^2/2 + 8*pi*a ./ x(i)));
        int_tf_func = @(k, tf) k.^2 .* E_k_f(k) ./ (exp(E_k_f(k) ./ tf) - 1);
        int_tf_func_term = @(tf) integral(@(k) int_tf_func(k, tf), 0, k_max, 'ArrayValued', true);
        
        % 定义方程
        equ = @(tf) int_ti_func_term / x(i) + 4 * pi^3 * a / x(i) * ...
            ((1 - 1 / x(i)) + 128 / (15 * sqrt(pi)) * a^(3/2) * ...
            (1 - (1 / x(i))^(3/2))) - int_tf_func_term(tf);
        
        % 求解 tf(i)
        options = optimset('Display', 'off', 'TolX', 1e-6);
        tf_results(j, i) = fzero(equ, [0, 10], options);
    end
end

% 归一化最终温度
tf_results = tf_results / T_c_i;

% 绘图
figure;
hold on;
colors = lines(length(ti_values)); % 使用不同颜色表示不同 ti

for j = 1:length(ti_values)
    plot(x, tf_results(j, :), 'LineWidth', 1, 'Color', colors(j, :));
end

xlabel('$V_f/V_i$', 'Interpreter', 'latex');
ylabel('$T_f/T_{c,i}$', 'Interpreter', 'latex');
legend(arrayfun(@(t) sprintf('$T_i=%.1f T_{c,i}$', t), ti_values/T_c_i, 'UniformOutput', false), ...
    'Interpreter', 'latex', 'Location', 'best');
title('Bogoliubov Theory (Numerical)');
box on;
hold off;

%% Isoenergy Curves (Fixed Vf / Vi, Varying Ti)
clear; clc;

% 基础参数
x = 2.6; % 固定体积比 Vf / Vi
a = 0.041; % 初始密度，单位 1/n_i^(1/3)
T_c_i = 2 * pi / zeta(3/2)^(2/3); % 初始临界温度
ti_values = linspace(0.2, 0.5, 50) * T_c_i; % 初始温度从 0.2 到 0.5 倍 T_c
k_max = Inf; % 数值积分的上限 

% 初始化存储变量
tf_results = zeros(size(ti_values)); % 存储每个初始温度下的 tf 值

% 循环：对每个初始温度 ti
for j = 1:length(ti_values)
    ti = ti_values(j); % 当前初始温度
    
    % 初始温度积分
    E_k_i = @(k) sqrt(k.^2 / 2 .* (k.^2 / 2 + 8 * pi * a));
    int_ti_func = @(k) k.^2 .* E_k_i(k) ./ (exp(E_k_i(k) ./ ti) - 1);
    int_ti_func_term = integral(int_ti_func, 0, k_max, 'ArrayValued', true);
    
    % 终态温度积分函数定义
    E_k_f = @(k) sqrt(k.^2 / 2 .* (k.^2 / 2 + 8 * pi * a / x));
    int_tf_func = @(k, tf) k.^2 .* E_k_f(k) ./ (exp(E_k_f(k) ./ tf) - 1);
    int_tf_func_term = @(tf) integral(@(k) int_tf_func(k, tf), 0, k_max, 'ArrayValued', true);
    
    % 定义能量守恒方程
    equ = @(tf) int_ti_func_term / x + 4 * pi^3 * a / x * ...
        ((1 - 1 / x) + 128 / (15 * sqrt(pi)) * a^(3/2) * ...
        (1 - (1 / x)^(3/2))) - int_tf_func_term(tf);
    
    % 求解 tf
    options = optimset('Display', 'off', 'TolX', 1e-6);
    tf_results(j) = fzero(equ, [0, 10], options);
end

% 归一化最终温度
tf_results = tf_results ./ ti_values; % T_f / T_i

% 绘图
figure;
plot(ti_values / T_c_i, tf_results, 'LineWidth', 1);
xlabel('$T_i / T_{c,i}$', 'Interpreter', 'latex');
ylabel('$T_f / T_i$', 'Interpreter', 'latex');
title('Bogoliubov Theory');
grid on;
box on;