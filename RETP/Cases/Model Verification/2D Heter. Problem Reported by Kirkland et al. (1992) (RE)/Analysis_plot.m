clear
clc

%% Data and plot
load('results_2DRE.mat')
v = [-1 -400];
figure
contour((ds.re_pw(2:end-1,2:end-1,3)./9810),v);
title('Water Pressure Head T = 12.5 day')
view(0,270)
xlabel('X')
ylabel('Z')
