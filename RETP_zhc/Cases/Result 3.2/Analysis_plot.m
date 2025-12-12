clear
clc

load('All_results.mat')
load("confi_process.mat")

for ii = 1:14
    r(ii).data = [par_vg(idxList{ii},:) sum1(idxList{ii})];
end

load('Result_supply.mat')
r(15).data = result_f2;

s= [];
for jj = 1:15
    s = [s; r(jj).data];
end

%% Figures

figure
hold on

% 定义蓝色渐变（从浅蓝到深蓝）
colors = [
    0.6 0.8 1;    % 浅蓝
    0.4 0.6 0.9;  % 中蓝
    0.2 0.4 0.8;  % 蓝
    0.1 0.3 0.7;  % 深蓝
    0.0 0.2 0.6    % 更深蓝
];

% 创建图例标签
legend_labels = cell(1, 5);

for ii = 15:-1:11
    idx = ii - 10;  % 颜色索引
    
    % 使用filled参数创建实心点，指定颜色
    scatter(r(ii).data(:,6)./r(ii).data(:,5), r(ii).data(:,7), 50, ...
            colors(idx, :), 'filled')
    
    set(gca, 'XScale', 'log')
    xlabel('k_{var}/k_{ref}')
    ylabel('RID')
    ylim([1e-4 1])
    
    % 创建图例标签
    ll = [1.0 0.8 0.6 0.4 0.2]';
    legend_labels{idx} = sprintf('R_H = %.1f', ll(ii-10));
    
    grid on
end

% 设置多行标题
title({'2D projection showing the relationship between k_{var}/k_{ref} and RID', ...
       'at T_{tot} under five spatial configurations when R_V = 1.0.'})

% 添加图例
legend(legend_labels, 'Location', 'best')

hold off

%% Rank Analysis

% 循环遍历 r(1) 到 r(15)
SRRC_weights = [];
R2_values = [];
for i = 1:15
    % 调用 rank_importance_two_soils_v2 函数并获取输出
    out = rank_importance_two_soils_v2(r(i).data);
    
    % 获取 SRRC_weight 和 SRRC_R2，并将它们添加到矩阵中
    SRRC_weights = [SRRC_weights; out.SRRC_weight];  % 每次迭代都将结果追加到矩阵中
    R2_values = [R2_values; out.SRRC_R2];
end

% 创建一个表格来显示所有结果
result_table = table(SRRC_weights, R2_values, 'VariableNames', {'SRRC_Weight', 'SRRC_R2'});

% 显示表格
disp(result_table);

data_matrix = table2array(result_table);

% 计算每列的统计量（均值、标准差、标准误）
means = mean(data_matrix, 1, 'omitnan');      % 每列的均值
stds = std(data_matrix, 0, 1, 'omitnan');     % 每列的标准差
n = sum(~isnan(data_matrix), 1);              % 每列的有效样本数

% 创建图形
figure('Position', [100, 100, 800, 600])
hold on

% 设置x轴位置
x = 1:size(data_matrix, 2);
bar_width = 0.7;

% 绘制柱状图
bar_handle = bar(x, means, bar_width, 'FaceColor', [0.2, 0.4, 0.8], 'EdgeColor', 'k', 'LineWidth', 1.5);

% 添加误差线（这里使用标准差作为误差线）
errorbar(x, means, stds, 'k.', 'LineWidth', 1.5, 'CapSize', 15);

% 或者使用标准误作为误差线（注释掉上面一行，使用下面一行）
% errorbar(x, means, stderr, 'k.', 'LineWidth', 1.5, 'CapSize', 15);

% 设置x轴标签
xtick_labels = {'$\alpha_{\mathrm{var}}/\alpha_{\mathrm{ref}}$', ...
                '$n_{\mathrm{var}}/n_{\mathrm{ref}}$', ...
                '$k_{\mathrm{var}}/k_{\mathrm{ref}}$', ...
                '$R^2$'};
set(gca, 'XTick', x, ...
         'XTickLabel', xtick_labels, ...
         'TickLabelInterpreter', 'latex', ...
         'FontSize', 14)
% 添加标签
xlabel('', 'FontSize', 14, 'FontWeight', 'bold')
ylabel('数值', 'FontSize', 14, 'FontWeight', 'bold')
title({'Standardized rank regression coefficient (SRRC^{2}) ', ...
    'and the corresponding R^{2} values.'}, 'FontSize', 16, 'FontWeight', 'bold')

% 添加网格
grid on
box on

hold off
