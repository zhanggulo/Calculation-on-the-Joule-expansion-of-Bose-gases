clear; clc; close all;

% 参数定义
a_values = linspace(0, 0.05, 100); % 密度范围 [0, 0.05]
k_max_1 = 10; % 截断积分上限，避免发散
k_max_2 = Inf; % 无穷积分上限
T_c = 2 * pi / zeta(3/2)^(2/3); % 临界温度
t = 0.1; % 固定温度 (单位：T_c)

% 色散关系
E_k = @(k, a) sqrt((k.^2 / 2) .* (k.^2 / 2 + 8 * pi * a));

% 分布函数
f = @(k, a, t) 1 ./ (exp(E_k(k, a) ./ t) - 1);

% 被积函数
int_func = @(k, a, t) k.^2 .* E_k(k, a) .* f(k, a, t);

% 数值积分项（有限上限）
int_term_1 = @(t, a) arrayfun(@(a) integral(@(k) int_func(k, a, t), 0, k_max_1, ...
    'AbsTol', 1e-6, 'RelTol', 1e-4), a);

% 数值积分项（无穷上限）
int_term_2 = @(t, a) arrayfun(@(a) integral(@(k) int_func(k, a, t), 0, k_max_2, ...
    'AbsTol', 1e-6, 'RelTol', 1e-4), a);

% 计算积分
int_val_1 = int_term_1(t * T_c, a_values);
int_val_2 = int_term_2(t * T_c, a_values);
error = abs(int_val_1 - int_val_2);

% 绘图
figure;
plot(a_values, error, 'LineWidth', 1);
xlabel('$n^{1/3}a$', 'Interpreter', 'latex');
ylabel('Integral Error');
grid on; box on;