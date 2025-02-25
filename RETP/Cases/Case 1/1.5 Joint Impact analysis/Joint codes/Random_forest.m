% Random Forest Out-of-Bag Permutation Error Analysis
clear
clc

% Input
load('All_results_mass.mat')
km = km * 0.01 / 86400 * 0.001 / 9.81 / 1000;
Rm = Rm * 0.01;
X = [nm am km];
Y = [Rm];

% Normalization
Xn = (X - mean(X)) ./ std(X);

%% Ten Repetition Analysis
imp_rec = [];
mse_rec = [];
r2_rec = [];

for ii = 1:10

% Set random forest parameters
numTrees = 500;  
minleaf = 5;

% Basic configuration (satisfies 85% of regular analysis needs)
model = TreeBagger(numTrees, Xn, Y, ...    % 500 trees balance efficiency and stability
    'Method', 'regression', ...
    'OOBPredictorImportance', 'on', ...
    'MinLeafSize', minleaf, ...           % Default value to avoid overfitting
    'PredictorSelection', 'curvature'); % Capture interaction effects

imp = model.OOBPermutedVarDeltaError;  % Feature importance
imp_rec = [imp_rec; imp];

% Validation
predictedY = predict(model, Xn);
% Compute Mean Squared Error (MSE)
mse = mean((Y - predictedY).^2);
mse_rec = [mse_rec; mse];

% Compute R² (Coefficient of Determination)
ssTotal = sum((Y - mean(Y)).^2);  % Total Sum of Squares
ssResidual = sum((Y - predictedY).^2);  % Residual Sum of Squares
r2 = 1 - (ssResidual / ssTotal);  % Compute R²
r2_rec = [r2_rec; r2];

end

% Take averages
m_r2 = mean(r2_rec);
m_imp = mean(imp_rec);

% Display results
disp(['The importance of n, α and k:', num2str(m_imp)]);
disp(['R²: ', num2str(r2)]);
