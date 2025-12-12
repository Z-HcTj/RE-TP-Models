clear
clc

%% Data
Output_name = {'1_Base_results.mat', '2_IC_Moist_results.mat', '3_IC_Wet_results.mat'};
% RE water content
re_thew_plot = [];
% TP water content
tp_thew_plot = [];
% RID
RID_index_plot = [];
% Mass balance 
re_wmass = [];
tp_wmass = [];

for ii = 1 : length(Output_name)
    load(Output_name{ii});
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

%% Water content profiles
z = [0.005:-0.01:-1.005]';
colors = [1 0.5 0; 0.53 0.81 0.98; 0 0.45 0.74];
re_wc_name = {"Base RE" "IC Moist RE" "IC Wet RE"};
tp_wc_name = {"Base TP" "IC Moist TP" "IC Wet TP"};
figure
hold on
for ii = 1:3
    plot(z,re_thew_plot(:,ii,end),'Color',colors(ii,:),'LineStyle','-','DisplayName',re_wc_name{ii})
    plot(z,tp_thew_plot(:,ii,end),'Color',colors(ii,:),'LineStyle','--','DisplayName',tp_wc_name{ii})
end
hold off
view(90,-90)
ylim([0 0.45])
xlim([-1 0])
legend
title('Water Content Profiles')
xlabel('Z [m]')
ylabel('Water Content')

%% RID plot
RID_name = {"Base" "IC Moist" "IC Wet"};
figure
hold on
for ii = 1:3
    plot(ds.trec,RID_index_plot(:,ii),'Color',colors(ii,:),'LineStyle','-','DisplayName',RID_name{ii})
end
hold off
legend
xlabel('Normalized Time')
ylabel('RID [%]')

