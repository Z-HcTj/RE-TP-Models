% Statistical analysis of data from Gupta et al. (2022), 
% plotting the distribution of MG parameters α and n, and the probability density of k.
clear
clc

%% Distribution of n and α
% Extract data from EXCEL
efname = 'Cleaned_WRC_Data.xlsx';
sheet = 1;
data_a = xlsread(efname, sheet, 'D:D');
data_n = xlsread(efname, sheet, 'E:E');
data_k = xlsread(efname, sheet, 'C:C');
data_n_a = [data_n data_a];
unique_data_n_a = unique(data_n_a, 'rows');

% Obtain the maximum α value corresponding to each n value and fit the curve
% Calculate the minimum interval between each point's n values
length_max = max(unique_data_n_a(2:end,1)-unique_data_n_a(1:end-1,1));
% Determine the number of separable n intervals
length_part = ceil((unique_data_n_a(end,1)-unique_data_n_a(1,1)) / 0.1);
% Get the maximum α value in each n interval
n_part_rangelef = unique_data_n_a(1,1);
n_part_rangerig = unique_data_n_a(1,1);
alpha_max_mat = []; % Array of maximum α values in each interval
n_max_mat = []; % Corresponding n values for maximum α
alpha_2max_mat = [];
n_2max_mat = [];
for ii = 1:length_part
    n_part_rangerig = n_part_rangerig + 0.1;
    part_range = unique_data_n_a(:,1) < n_part_rangerig & unique_data_n_a(:,1) >= n_part_rangelef;
    sorted_a_values = sort(unique_data_n_a(part_range,2), 'descend');
    alpha_max = sorted_a_values(1);
    n_max = unique_data_n_a(find(unique_data_n_a(:,2) == alpha_max),1);
    alpha_max_mat = [alpha_max_mat; alpha_max];
    n_max_mat = [n_max_mat; n_max];
    if length(unique_data_n_a(part_range,2)) > 1
        alpha_2max = sorted_a_values(2);
        n_2max = unique_data_n_a(find(unique_data_n_a(:,2) == alpha_2max),1);
        alpha_2max_mat = [alpha_2max_mat; alpha_2max];
        n_2max_mat = [n_2max_mat; n_2max];
    end
    n_part_rangelef = n_part_rangerig;
end
% Combine the results
alpha_max_mat = [alpha_max_mat; alpha_2max_mat];
n_max_mat = [n_max_mat; n_2max_mat];
% Optimize the data corresponding to n values
n_ft_range = find(n_max_mat>1.5);
n_ft_mat = n_max_mat(n_ft_range);
alpha_ft_mat = alpha_max_mat(n_ft_range);
% Starting point
start_point = [1.5, 100];
% Fit parameters
x = n_ft_mat; 
y = alpha_ft_mat; 
% Weight vector, assign higher weights to areas with fewer data points
weights = ones(size(x));
weights(x<1.6) = 1;
% Define custom fitting function
ft = fittype('a*x^(-b) + 4', 'independent', 'x', 'coefficients', {'a', 'b'});
% Fit the data using custom weights
[fitresult, gof] = fit(x, y, ft, 'Weights', weights);

% Plot the distribution of MG parameters α and n, along with the fitted curve
figure
scatter(data_n_a(:,1),data_n_a(:,2))
grid on
hold on
scatter(n_max_mat,alpha_max_mat,10,'r')
scatter(n_max_mat(10:end),alpha_max_mat(10:end),10,'r')
xFit = [1.5:0.01:7]';
yFit = fitresult(xFit);
plot(xFit, yFit, 'r-', 'LineWidth', 2);
hold off;
xlabel("n")
ylabel("α")

% Output the fitting parameters
fprintf('Fitting parameter a = %f\n', fitresult.a);
fprintf('Fitting parameter b = %f\n', fitresult.b);
% Goodness of fit information
disp(gof)

%% Plot the probability density of k
% Probability distribution of k
indices = find(data_n > 1.5);
selected_n = data_n(indices);
selected_a = data_a(indices);
selected_k = data_k(indices);
binEdges = 0:1:8000;
[counts, edges] = histcounts(selected_k, binEdges);
binCenters = edges(1:end-1) + diff(edges) / 2;
binCenters = binCenters';
counts = counts' / 5251;
% Calculate Q-90 and Q-10
q90 = prctile(selected_k, 90); % 90th percentile
q10 = prctile(selected_k, 10); % 10th percentile
% Plotting
figure;
histogram(selected_k, 'BinEdges', binEdges, 'Normalization', 'pdf', 'EdgeColor', 'none'); 
xlabel('Ks values (0-100 cm/d)');
ylabel('Probability Density');
xlim([0 100]); % Limit the x-axis range
grid on;
