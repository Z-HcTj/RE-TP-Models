clear
clc

%% Data
load("Case_name.mat")
op = 0;
% RE water content
re_thew_plot = [];
% TP water content
tp_thew_plot = [];
% RID
RID_index_plot = [];
re0 = [];
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
    re0 = [re0; (re_sum(end) - re_sum(1))];
    RID_index_plot = [RID_index_plot RID_index];
    re_wmass = [re_wmass ds.re_wmass];
    tp_wmass = [tp_wmass ds.tp_wmass];
end
RID_index_plot = RID_index_plot * 100;

colors = [0 0.45 0.74; 0, 0.7450, 0.9330; 0.53 0.81 0.98; 0.98 0.85 0.73; 1 0.5 0; 1 0.27 0];
RID_name = Output_name;
figure
hold on
for ii = 1:6
    plot(ds.trec,RID_index_plot(:,ii),'Color',colors(ii,:),'LineStyle','-','DisplayName',RID_name{ii})
end
hold off
legend
xlabel('Normalized Time')
ylabel('RID [%]')

%% uRID
% RE water content
re_thew_plot = [];
% TP water content
tp_thew_plot = [];
% RID
uRID_index_plot = [];
tt = 42;

for ii = 1 : length(Output_name)
    load(Output_name{ii});
    re_thew_plot = [re_thew_plot ds.re_thew(:,2,:)];
    tp_thew_plot = [tp_thew_plot ds.tp_thew(:,2,:)];
    re_sum = [];
    tp_sum = [];
    for jj = 1 : length(ds.trec)
    re_seff = (re_thew_plot(2:tt-1,ii,jj) - ds.thewr(2:tt-1,2)) ...
        ./ (ds.phi(2:tt-1,2) - ds.thewr(2:tt-1,2));
    tp_seff = (tp_thew_plot(2:tt-1,ii,jj) - ds.thewr(2:tt-1,2)) ...
        ./ (ds.phi(2:tt-1,2) - ds.thewr(2:tt-1,2));
    re_sum = [re_sum; sum(re_seff)];
    tp_sum = [tp_sum; sum(tp_seff)];
    end
    RID_index = (re_sum - tp_sum) ./ (re_sum(end) - re_sum(1));
    if op == 1
    RID_index = (re_sum - tp_sum) ./ re0(ii);
    end
    uRID_index_plot = [uRID_index_plot RID_index];
    re_wmass = [re_wmass ds.re_wmass];
    tp_wmass = [tp_wmass ds.tp_wmass];
end
uRID_index_plot = uRID_index_plot * 100;

colors = [0 0.45 0.74; 0, 0.7450, 0.9330; 0.53 0.81 0.98; 0.98 0.85 0.73; 1 0.5 0; 1 0.27 0];
RID_name = Output_name;
figure
hold on
for ii = 1:6
    plot(ds.trec,uRID_index_plot(:,ii),'Color',colors(ii,:),'LineStyle','-','DisplayName',RID_name{ii})
end
hold off
legend
xlabel('Normalized Time')
ylabel('uRID [%]')

%% dRID plot
% RE water content
re_thew_plot = [];
% TP water content
tp_thew_plot = [];
% RID
dRID_index_plot = [];
tt = 42;

for ii = 1 : length(Output_name)
    load(Output_name{ii});
    re_thew_plot = [re_thew_plot ds.re_thew(:,2,:)];
    tp_thew_plot = [tp_thew_plot ds.tp_thew(:,2,:)];
    re_sum = [];
    tp_sum = [];
    for jj = 1 : length(ds.trec)
    re_seff = (re_thew_plot(tt:end-1,ii,jj) - ds.thewr(tt:end-1,2)) ...
        ./ (ds.phi(tt:end-1,2) - ds.thewr(tt:end-1,2));
    tp_seff = (tp_thew_plot(tt:end-1,ii,jj) - ds.thewr(tt:end-1,2)) ...
        ./ (ds.phi(tt:end-1,2) - ds.thewr(tt:end-1,2));
    re_sum = [re_sum; sum(re_seff)];
    tp_sum = [tp_sum; sum(tp_seff)];
    end
    RID_index = (re_sum - tp_sum) ./ (re_sum(end) - re_sum(1));
    if op == 1
    RID_index = (re_sum - tp_sum) ./ re0(ii);
    end
    dRID_index_plot = [dRID_index_plot RID_index];
    re_wmass = [re_wmass ds.re_wmass];
    tp_wmass = [tp_wmass ds.tp_wmass];
end
dRID_index_plot = dRID_index_plot * 100;

colors = [0 0.45 0.74; 0, 0.7450, 0.9330; 0.53 0.81 0.98; 0.98 0.85 0.73; 1 0.5 0; 1 0.27 0];
RID_name = Output_name;
figure
hold on
for ii = 1:6
    plot(ds.trec,dRID_index_plot(:,ii),'Color',colors(ii,:),'LineStyle','-','DisplayName',RID_name{ii})
end
hold off
legend
xlabel('Normalized Time')
ylabel('dRID [%]')

