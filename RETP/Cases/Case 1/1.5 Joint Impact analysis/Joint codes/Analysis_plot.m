clear
clc

%% Data
load('All_results.mat')
load('LHS_points.mat')
result_suffix = 1:1000;
% RE water content
re_thew_plot = [];
% TP water content
tp_thew_plot = [];
% RID
RID_index_plot = [];
% Mass balance 
re_wmass = [];
tp_wmass = [];

for ii = 1 : length(result_suffix)
    ds = All_results{1,ii}.ds;
    re_thew_plot = [re_thew_plot ds.re_thew(:,2,:)];
    tp_thew_plot = [tp_thew_plot ds.tp_thew(:,2,:)];
    re_sum = [];
    tp_sum = [];
    for jj = 1 : length(ds.trec)
    re_seff = (re_thew_plot(2:end-1,ii,jj) - ds.thewr(2:end-1,2)) ...
        ./ (ds.phi(2:end-1,2) - ds.thewr(2:end-1,2));
    tp_seff = (tp_thew_plot(2:end-1,ii,jj) - ds.thewr(2:end-1,2)) ...
        ./ (ds.phi(2:end-1,2) - ds.thewr(2:end-1,2));
    re_sum = [re_sum; sum(re_seff)];
    tp_sum = [tp_sum; sum(tp_seff)];
    end
    RID_index = (re_sum - tp_sum) ./ (re_sum(end) - re_sum(1));
    RID_index_plot = [RID_index_plot RID_index];
    re_wmass = [re_wmass ds.re_wmass];
    tp_wmass = [tp_wmass ds.tp_wmass];
end
RID_index_plot = RID_index_plot * 100;

%% RID plot
% Remove samples that do not satisfy mass conservation
re_mass_find = find(re_wmass > 0.001);
tp_mass_find = find(tp_wmass > 0.001);
mass_find = [re_mass_find tp_mass_find];
mass_find = unique(mass_find);
mass_right = setdiff(1:length(All_results), mass_find);

% Plot
figure
scatter3(sampled_params.n(mass_right),sampled_params.alpha(mass_right),RID_index_plot(end,mass_right))
xlabel('n')
ylabel('α [1/m]')
zlabel('RID [%]')

% Save for ''Random Forest''
nm = sampled_params.n(mass_right);
am = sampled_params.alpha(mass_right);
km = sampled_params.k(mass_right);
Rm = [RID_index_plot(end,mass_right)]';
save("All_results_mass.mat", 'nm', 'am'...
    ,'km', 'Rm')

