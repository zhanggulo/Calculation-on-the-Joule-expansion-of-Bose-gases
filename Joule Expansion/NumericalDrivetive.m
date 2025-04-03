clear; clc;

% 参数设置
n = 0.03 ^ 3;          % 固定密度
t = 0.1;           % 固定温度
delta_n = 1e-9;    % 密度微小变化量
k = linspace(1e-6, 10, 1000); % k 的范围，避免 k=0 处的数值问题

% 色散关系
E_k = @(k, n) sqrt((k.^2 / 2) .* (k.^2 / 2 + 8 * pi * n));

% int_func_1
int_func_1 = @(k, n, t) 4 * pi .* (k.^2 ./ 2) ./ E_k(k, n) .* ...
    (exp(E_k(k, n)) .* (1 - E_k(k, n)./t) - 1) ./ (exp(E_k(k, n) ./ t) - 1).^2;

% 分布函数
f = @(k, n, t) E_k(k, n) ./ (exp(E_k(k, n) ./ t) - 1);

% df_dn (数值微分)
df_dn = @(k, n, t) (f(k, n + delta_n, t) - f(k, n - delta_n, t)) / (2 * delta_n);

% 计算 int_func_1 和 df_dn
int_func_1_vals = arrayfun(@(k_val) int_func_1(k_val, n, t), k);
df_dn_vals = arrayfun(@(k_val) k_val.^2 .* df_dn(k_val, n, t), k);

error = abs(int_func_1_vals - df_dn_vals);

% 绘图
figure;
plot(k, error, 'LineWidth', 1);
xlabel('Wavevector k');
ylabel('Function Values');
title('Comparison of int\_func\_1 and df\_dn vs k');
grid on;

clear; clc; close all;

% 参数设置
n_values = linspace(0, 0.05, 50).^3; % 密度范围 [0.001, 0.05]^3
t = 0.1;            % 固定温度
delta_n = 1e-9;     % 密度微小变化量
k = linspace(1e-6, 10, 1000); % 避免 k=0 的数值问题

% 色散关系
E_k = @(k, n) sqrt((k.^2 / 2) .* (k.^2 / 2 + 8 * pi * n));

% int_func_1 (解析积分)
int_func_1 = @(k, n, t) 4 * pi .* (k.^2 ./ 2) ./ E_k(k, n) .* ...
    (exp(E_k(k, n)./t) .* (1 - E_k(k, n)./t) - 1) ./ (exp(E_k(k, n) ./ t) - 1).^2;

% 分布函数
f = @(k, n, t) E_k(k, n) ./ (exp(E_k(k, n) ./ t) - 1);

% 数值微分 df/dn
df_dn = @(k, n, t) (f(k, n + delta_n, t) - f(k, n - delta_n, t)) / (2 * delta_n);

% 初始化误差数组
integral_error = zeros(size(n_values));

% 循环计算不同 n 下的积分误差
for i = 1:length(n_values)
    n = n_values(i);
    
    % 计算解析积分
    int_func_1_vals = arrayfun(@(k_val) k_val.^2 .* int_func_1(k_val, n, t), k);
    int_analytical = trapz(k, int_func_1_vals);
    
    % 计算数值积分
    df_dn_vals = arrayfun(@(k_val) k_val.^2 .* df_dn(k_val, n, t), k);
    int_numerical = trapz(k, df_dn_vals);
    
    % 计算误差
    integral_error(i) = abs(int_analytical - int_numerical);
end

% 绘制误差图
figure;
plot(n_values.^(1/3), integral_error, '-r', 'LineWidth', 1.5);
xlabel('$n^{1/3}$', 'Interpreter', 'latex');
ylabel('Absolute Integral Error');
title('Integral Error as a Function of $n^{1/3}$');
grid on; box on;