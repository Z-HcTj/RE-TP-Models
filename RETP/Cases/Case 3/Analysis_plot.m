clear
clc

%% Data
Output_name = {'2DHeter.mat'};
load('2DHeter.mat')
% RE seff
re_seff_plot = (ds.re_thew - ds.thewr) ./ (ds.phi - ds.thewr);
% TP seff
tp_seff_plot = (ds.tp_thew - ds.thewr) ./ (ds.phi - ds.thewr);
% RE-TP seff
re_tp_seff_plot = re_seff_plot - tp_seff_plot;

% T = 29000s
figure
heatmap(re_tp_seff_plot(2:end-1,2:end-1,1))
ax = gca;
ax.XDisplayLabels = nan(size(ax.XDisplayData));
ax.YDisplayLabels = nan(size(ax.YDisplayData));
title('T = 29000 s')

% T = 58000s
figure
heatmap(re_tp_seff_plot(2:end-1,2:end-1,2))
ax = gca;
ax.XDisplayLabels = nan(size(ax.XDisplayData));
ax.YDisplayLabels = nan(size(ax.YDisplayData));
title('T = 58000 s')

