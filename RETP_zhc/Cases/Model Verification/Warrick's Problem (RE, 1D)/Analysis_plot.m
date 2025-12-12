clear
clc

%% Data
Output_name = {'results_Warrick.mat'};
% RE water content
re_thew_plot = [];

for ii = 1 : length(Output_name)
    load(Output_name{ii});
    re_thew_plot = [re_thew_plot ds.re_thew(:,2,:)];
end

%% Water content profiles
z = [0.005:-0.01:-1.005]';
colors = ['b' 'r' 'k'];
re_wc_name = {"RE model"};
% Analytical solution
wc1 = [0.0825 0.0825 0.0825 ;0.165 0.165 0.165 ;0.2475 0.2475 0.2475];
z1 = [74.58 61.32 39.02; 75.20 62.13 40.04; 77.12 64.54 42.80] - 100;
z1 = z1 * 0.01;

figure
hold on
for ii = 1:3
    plot(z,re_thew_plot(:,:,ii),'Color',colors(ii),'LineStyle','-')
    plot(z1(ii,:),wc1(ii,:),'*','Color',colors(ii))
end
hold off
view(90,-90)
xlim([-1 0])
title('Water Content Profiles')
xlabel('Z [m]')
ylabel('Water Content')
