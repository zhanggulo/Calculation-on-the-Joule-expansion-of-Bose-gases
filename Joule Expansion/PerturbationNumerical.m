clc; clear;

% 参数
t_values = 0.15;  % 初始无量纲温度数组
n1_3_a = 0.02;
v_values = linspace(1, 8, 1000);  % 细化 v 取值范围

% 计算前系数
C = (2 * zeta(3/2)^(5/3)) / (3 * zeta(5/2)) * n1_3_a;  

% figure;
hold on;

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
        t_f_initial = max(t, 0.01);  % 避免 t_f 过小，导致数值不稳定

        % 求解 t_f
        t_f_values(i) = fsolve(fun, t_f_initial, optimoptions('fsolve','Display','off')) / t;
    end

    % 绘制曲线
    plot(v_values, t_f_values, 'b-', 'LineWidth', 0.6);

    % 计算极大值点
    [t_f_max, max_index] = max(t_f_values);
    v_max = v_values(max_index);
    plot(v_max, t_f_max, 'r^', 'MarkerSize', 6, 'LineWidth', 0.6); % 红色三角形标记极大值

    % 找到最接近 1 的点（但排除 v = 1 附近）
    valid_indices = find(v_values > 1.1);  % 只考虑 v > 1.1
    if ~isempty(valid_indices)
        [~, one_index] = min(abs(t_f_values(valid_indices) - 1));
        one_index = valid_indices(one_index); % 获取原始索引
        plot(v_values(one_index), t_f_values(one_index), 'ro', 'MarkerSize', 6, 'LineWidth', 0.6); % 红色圆圈
    end
end

% 图像设置
xlabel('$V_f / V$', 'Interpreter', 'latex', 'FontSize', 10);
ylabel('$T_f/T$','Interpreter', 'latex', 'FontSize', 10);
hold off;
box on;

%%
% 参数
T_values = linspace(0.1, 0.5, 1000);  % 无量纲温度数组
n1_3_a = 0.02;
v_values = 5;  % 固定的 v 取值

% 计算前系数
C = (2 * zeta(3/2)^(5/3)) / (3 * zeta(5/2)) * n1_3_a;  

% figure;
hold on;

% 遍历不同的 v 进行求解
for j = 1:length(v_values)
    v = v_values(j);  % 取当前 v
    t_f_values = zeros(size(T_values)); % 存储解

    % 遍历不同的 T 进行数值求解
    for i = 1:length(T_values)
        T = T_values(i);

        % 定义非线性方程
        fun = @(t_f) t_f^(5/2) + C * (2*t_f^3 - t_f^(3/2) / v) - ...
                     (T^(5/2) / v + C * (1 - 1/v - T^(3/2) + 2*T^3) / v);

        % 设定初始猜测值 t_f_initial
        t_f_initial = max(T, 0.01);  % 避免 t_f 过小，导致数值不稳定

        % 求解 t_f
        t_f_values(i) = fsolve(fun, t_f_initial, optimoptions('fsolve','Display','off')) / T;
    end

    % 绘制曲线
    plot(T_values, t_f_values, 'b-', 'LineWidth', 0.6);

    % 找到最接近 1 的点并标记
    [~, one_index] = min(abs(t_f_values - 1)); % 找到误差最小的点
    plot(T_values(one_index), t_f_values(one_index), 'ro', 'MarkerSize', 6, 'LineWidth', 0.6); % 红色圆圈
end

% 图像设置
xlabel('$T/T_c$', 'Interpreter', 'latex', 'FontSize', 10);
ylabel('$T_f/T$','Interpreter', 'latex', 'FontSize', 10);
hold off;
box on;
