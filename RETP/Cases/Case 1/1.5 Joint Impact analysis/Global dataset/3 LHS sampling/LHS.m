% Use the LHS (Latin Hypercube Sampling) method to randomly sample n, α, and k values within a specified range.
% The results of each sampling will differ every time the program is run.
clear
clc

%% LHS Sampling Points
% Total number of samples
num_samples = 1000;

% Latin Hypercube Sampling
lhs_samples = lhsdesign(num_samples, 3);

% Range of n
n_min = 1.5;
n_max = 7;
sampled_n = n_min + (n_max - n_min) * lhs_samples(:, 1);

% Dynamically adjust the range of α
alpha_max_func = @(n) 177.640847 * n.^(-3.641711) + 4;
sampled_alpha = arrayfun(@(n_val, alpha_norm) ...
    max(0.5, alpha_max_func(n_val) * alpha_norm), ...
    sampled_n, lhs_samples(:, 2));

% Range of k
k_min = 1;
k_max = 400;
sampled_k = k_min + (k_max - k_min) * lhs_samples(:, 3);

% 3D Scatter Plot
figure
scatter3(sampled_n, sampled_alpha, sampled_k, 36, sampled_k, 'filled')
xlim([1.5, 7])
ylim([0.5, 40])
zlim([1, 400])
colorbar
title('3D Sampling Distribution')
xlabel('n')
ylabel('alpha')
zlabel('K')

% Generate a table of sampled results
sampled_params = table(sampled_n, sampled_alpha, sampled_k,...
    'VariableNames', {'n', 'alpha', 'k'});
