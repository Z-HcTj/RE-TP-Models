% =========================================================================
% Check_Correlation_Fixed.m
% 用于计算 n, alpha, k 之间的皮尔逊相关系数 (含数据类型自动修复)
% =========================================================================

clc; clear;

filename = 'LHS_points.mat'; 

if exist(filename, 'file')
    load(filename);
    
    % 检查变量是否存在
    if ~exist('sampled_params', 'var')
        error('文件中没找到变量 sampled_params，请检查变量名！');
    end
    
    raw_data = sampled_params;
    
    % --- 关键修复步骤：数据类型转换 ---
    disp(['检测到数据类型: ', class(raw_data)]);
    
    if istable(raw_data)
        % 如果是 Table，转换为 Array
        disp('正在将 Table 转换为数值矩阵...');
        data = table2array(raw_data);
    elseif isstruct(raw_data)
        % 如果是 Struct，尝试提取字段 (需要根据你的实际字段名修改)
        error('数据是结构体 (struct)，请手动提取三列数据，例如: data = [raw_data.n, raw_data.alpha, raw_data.k];');
    elseif iscell(raw_data)
        % 如果是 Cell，转换为 Matrix
        disp('正在将 Cell 转换为数值矩阵...');
        data = cell2mat(raw_data);
    else
        % 默认假设已经是矩阵
        data = raw_data;
    end
    
    % 再次检查是否为纯数值
    if ~isnumeric(data)
        error('转换失败：数据包含非数值内容，请检查源文件。');
    end

    % 检查维度
    [rows, cols] = size(data);
    fprintf('数据加载成功：共 %d 行，%d 列\n', rows, cols);
    
    if cols ~= 3
        error('数据列数不对！必须是 3 列 (n, alpha, k)。检测到 %d 列。', cols);
    end
else
    error('未找到文件 LHS_points.mat，请确保文件在当前路径下。');
end

% 2. 计算相关系数矩阵
R = corrcoef(data);

% 提取三个两两配对的相关系数
% 假设列顺序是: 1->n, 2->alpha, 3->k
r_n_alpha = R(1, 2);
r_n_k     = R(1, 3);
r_alpha_k = R(2, 3);

% 3. 输出结果到屏幕
fprintf('\n--------------------------------------------------\n');
fprintf('皮尔逊相关系数 (Pearson Correlation Coefficients):\n');
fprintf('--------------------------------------------------\n');
fprintf('n 与 alpha 的相关系数 (r): %.4f\n', r_n_alpha);
fprintf('n 与 k     的相关系数 (r): %.4f\n', r_n_k);
fprintf('alpha 与 k 的相关系数 (r): %.4f\n', r_alpha_k);
fprintf('--------------------------------------------------\n');

% 4. 计算绝对值的最大值
max_abs_r = max(abs([r_n_alpha, r_n_k, r_alpha_k]));

fprintf('【论文填写指南】\n');
fprintf('正文中建议填写的最大相关性数值 (|r| < X): %.2f\n', max_abs_r);

if max_abs_r < 0.05
    fprintf('结论: 完美！参数高度独立 (Independent)。\n');
elseif max_abs_r < 0.3
    fprintf('结论: 很好！存在极弱相关 (Weak correlation)，是物理约束导致的正常现象。\n');
else
    fprintf('结论: 注意！存在较强相关，必须在回复信中强调这是物理约束 (alpha_max) 导致的。\n');
end