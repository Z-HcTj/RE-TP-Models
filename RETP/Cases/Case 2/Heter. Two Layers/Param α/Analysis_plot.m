clear
clc

%% Data
load("Case_name.mat")
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

%% RID plot
colors = [0 0.45 0.74; 0.53 0.81 0.98; 0.98 0.85 0.73; 1 0.5 0; 1 0.27 0];
RID_name = Output_name;
figure
hold on
for ii = 1:5
    plot(ds.trec,RID_index_plot(:,ii),'Color',colors(ii,:),'LineStyle','-','DisplayName',RID_name{ii})
end
hold off
legend
xlabel('Normalized Time')
ylabel('RID [%]')

delete("Case_name.mat")
