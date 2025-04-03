clc; clear;

% 参数
t_values = 0 : 0.01 : 0.3;  % 初始无量纲温度数组
n1_3_a = 0.02;
v_values = linspace(1, 8, 1000);

% 计算前系数
C = (2 * zeta(3/2)^(5/3)) / (3 * zeta(5/2)) * n1_3_a;  

% figure;
hold on;

% 存储峰值信息
t_f_matrix = zeros(length(t_values), length(v_values));
peak_v = zeros(size(t_values));
peak_t_f = zeros(size(t_values));

% 遍历不同的 t 进行求解
for j = 1:length(t_values)
    t = t_values(j);  % 取当前 t
    t_f_values = zeros(size(v_values)); % 存储解

    % 遍历不同的 v 进行数值求解
    for i = 1:length(v_values)
        v = v_values(i);

        % 定义非线性方程
        fun = @(t_f) t_f^(5/2) + C * (2*t_f^3 - t_f^(3/2) / v) - ...
                     (t^(5/2) / v + C * (1 - 1/v - t^(3/2) + 2*t^3) / v);

        % 设定初始猜测值 t_f_initial
        t_f_initial = 0.1;  % 避免 t_f 过小，导致数值不稳定

        % 求解 t_f
        t_f_values(i) = fsolve(fun, t_f_initial, optimoptions('fsolve','Display','off'));
    end

    % 存储数据
    t_f_matrix(j, :) = t_f_values;
    
    % 绘制曲线
    % plot(v_values, t_f_values, 'LineWidth', 0.6);
end

% 计算等能量曲线的最大值
for j = 1:length(t_values)
    [peak_t_f(j), peak_idx] = max(t_f_matrix(j, :));
    peak_v(j) = v_values(peak_idx);
end

% 绘制峰值曲线
% figure;
plot(peak_v, t_values, 'k--', 'LineWidth', 0.6);

% 图像设置
xlabel('$V_f / V$', 'Interpreter', 'latex', 'FontSize', 10);
ylabel('$T/T_c$','Interpreter', 'latex', 'FontSize', 10);
hold off;
box on;
