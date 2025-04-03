%% Isoenergy Curves (Ti/Tc was fixed)
clear; clc

eta = 1; % Condensation Fraction
x = linspace(1, 6, 1000); % Expansion fraction Vf / Vi
ni = 0.02; % Initial density in unit of (1/a)^3

% Perturbation 计算 t_i
t_i = (2 * zeta(3/2)^(5/3) / (3 * zeta(5/2)) * (2 - eta^2) * ni ./ x) .^ (2/5);

% 计算等能量曲线
ti = 0:0.01:0.3; % 初始温度范围 (单位 Tc)
tf = cell(size(ti));

for i = 1:length(ti)
    tf{i} = ((ti(i)) .^ (5/2) ./ x + 2 .* zeta(3/2) .^ (5/3) / (3 * zeta(5/2)) .* ...
        (2 - eta.^2) .* ni ./ x .* (1 - 1 ./ x)) .^ (2/5);
end

% 创建图像
figure
hold on

% 填充下右半部分（淡红色，透明）
fill([x, fliplr(x)], [t_i, fliplr(zeros(size(t_i)))], [1, 0.7, 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.8)

% 填充上半部分（淡蓝色，透明）到顶部
fill([x, fliplr(x)], [t_i, fliplr(ones(size(t_i)))], [0.6, 0.8, 1], 'EdgeColor', 'none', 'FaceAlpha', 0.8)

% 添加边界虚线
plot(x, t_i, 'k--', 'LineWidth', 0.6) % Perturbation 线（T=Ti的分界线）

% 重新绘制 T_i 主要曲线
plot(x, t_i, "k", "LineWidth", 0.6) 

% 坐标轴与标签
xlabel('$V_f/V_i$', 'Interpreter','latex')
ylabel('$T/T_{c,i}$','Interpreter','latex')
ylim([0 0.4])
hold off

%% 填充虚线右侧区域为透明的深红色
hold on;
fill([peak_v, fliplr(peak_v)], [t_values, fliplr(zeros(size(t_values)))], [0.95, 0.6, 0.6], 'EdgeColor', 'none', 'FaceAlpha', 0.8)

% 添加一个0.5宽的边框
set(gca, 'LineWidth', 0.5) % Set the box border line width
box on % Display the box around the plot
hold off;