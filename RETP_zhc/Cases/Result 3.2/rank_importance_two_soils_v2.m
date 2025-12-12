function out = rank_importance_two_soils_v2(data) 
% data 是一个 N x 7 的矩阵，包含基准土和相对土的α、n、k，以及模拟结果y
% 数据列数：
% 1-2: α（基准土、相对土）；3-4: n（基准土、相对土）；5-6: k（基准土、相对土）；7: 模拟结果y

% 提取各列数据
alpha1 = data(:, 1); alpha2 = data(:, 2); % 基准土与相对土的 α
n1 = data(:, 3); n2 = data(:, 4);         % 基准土与相对土的 n
k1 = data(:, 5); k2 = data(:, 6);         % 基准土与相对土的 k
y = data(:, 7);                           % 模拟结果 y

% ---- 1) 计算相对特征（log 比）----
rA = alpha2(:)./alpha1(:);  % 相对α差异
rN = n2(:)./n1(:);          % 相对n差异
rK = k2(:)./k1(:);          % 相对k差异

% ---- 2) 数据筛选（如果需要） ----
% 例如筛选出小于一定值的结果
% tt = find(rK < 1);  
% rA = rA(tt);
% rN = rN(tt);
% rK = rK(tt);
% y = y(tt);
% n1 = n1(tt);
% n2 = n2(tt);
% alpha1 = alpha1(tt);
% alpha2 = alpha2(tt);
% k1 = k1(tt);
% k2 = k2(tt);

% ---- 3) 将比值作为特征 ----
X = [rA, rN, rK];  % 只保留rA、rN、rK作为自变量矩阵
N = numel(y);  % 样本数量

% ---- 4) SRRC：秩变换 + 标准化回归 ----
Xr = zeros(N, 3);  % 3个特征（只包括rA, rN, rK）
for j = 1:3
    Xr(:, j) = tiedrank(X(:, j));  % 对特征进行秩变换
end
Yr = tiedrank(y(:));  % 模拟结果 y 的秩变换
Xz = zscore(Xr);      % 对 X 进行标准化
Yz = zscore(Yr);      % 对 y 进行标准化

% ---- 5) OLS 回归 ----
B = [ones(N, 1) Xz] \ Yz;  % OLS 回归
beta = B(2:end);           % 标准化系数
yhat = [ones(N, 1) Xz] * B;  % 拟合值
R2 = var(yhat) / var(Yz);  % 决定系数（拟合优度）

% ---- 6) 计算每个特征的贡献权重 ----
w = beta.^2;  % 计算权重（标准化系数的平方）
w = w / sum(w);  % 归一化为百分比

% ---- 7) Unique 贡献（残差法）----
Xpure = Xz;  % 3 列
uniqueR2 = zeros(1, 3);
for j = 1:3
    keep = setdiff(1:3, j);
    Bj = [ones(N, 1) Xpure(:, keep)] \ Xpure(:, j);  % 保留其中两个特征，计算残差
    rj = Xpure(:, j) - [ones(N, 1) Xpure(:, keep)] * Bj;  % 计算去相关后的“纯净 j”
    bj = [ones(N, 1) rj] \ Yz;  % 纯净部分单独回归
    yj = [ones(N, 1) rj] * bj;  % 预测值
    uniqueR2(j) = var(yj) / var(Yz);  % 计算唯一贡献的 R2
end
unique_share = uniqueR2 / sum(uniqueR2);  % 唯一贡献的归一化百分比

% ---- 8) 非参数校验（随机森林置换重要性）----
try
    M = TreeBagger(300, X, y(:), 'Method', 'regression', ...
        'OOBPrediction', 'On', 'OOBPredictorImportance', 'On', ...
        'PredictorNames', {'rAlpha', 'rN', 'rK'});
    permImp = M.OOBPermutedPredictorDeltaError(:)';  % 获取置换重要性
    permImp = permImp / sum(permImp);                % 归一化
catch
    permImp = nan(1, 3);  % 如果随机森林不可用，输出 NaN
end

% ---- 9) 自举置信区间（可选，稳健性）----
Bboot = 500;  % 自举次数
bootW = zeros(Bboot, 3);
rng(1);
for b = 1:Bboot
    idx = randi(N, N, 1);  % 随机抽样
    Xb = Xz(idx, :); Yb = Yz(idx);
    Bb = [ones(N, 1) Xb] \ Yb;
    betab = Bb(2:end);
    wb = betab.^2; wb = wb / sum(wb);  % 计算权重
    bootW(b, :) = wb;
end
ci = quantile(bootW, [0.025 0.975]);  % 计算置信区间

% ---- 10) 输出包装 ----
out.labels = {'rAlpha', 'rN', 'rK'};  % 输出特征标签
out.SRRC_weight = w(:)';  % 总贡献（秩回归
out.SRRC_R2 = R2;  % 决定系数
out.Unique_share = unique_share;  % 唯一贡献占比
out.Permutation = permImp;  % 非参数重要性（可能为 NaN）
out.SRRC_CI_95 = ci;  % 自举 95% CI
out.UR2 = uniqueR2;
end

