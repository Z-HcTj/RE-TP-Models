%% Run the program to get the output
clear
clc

% Path to ''main'' file
oper_filepath = fullfile(fileparts(fileparts(fileparts(fileparts(pwd)))), 'Main');
% Path back to this file
current_filepath = pwd;
save(fullfile(oper_filepath,'curpath.mat'),'current_filepath')

%% Single Impact
Output_name = {'Heter. k=0.236e-12.mat' 'Heter. k=0.708e-12.mat' 'Heter. k=1.77e-12.mat'...
    'Heter. k=4.67e-12.mat' 'Heter. k=9.64e-12.mat'};
key_param = [20 60 150 396 800];
save("Case_name","Output_name","key_param")

for ii = 1:length(key_param)
    save("run_ii","ii")

    % Input parameters for every scenario
    Inputdata

    % Save
    ds = struct();
    ds.info = 'Output Results';
    save(fullfile(oper_filepath,'mydata.mat'), 'ds');

    % To ''main'' file and run
    cd(oper_filepath)

    % Get total simulation time
    re_gettime

    % Get output results
    re_main
    tp_main

    % Output results transfer to this file
    load('curpath.mat')
    save(fullfile(current_filepath,'mydata.mat'), 'ds');
    cd(current_filepath)
    load('Case_name.mat')
    load('run_ii.mat')
    movefile('mydata.mat',  Output_name{ii}); 
    
end

delete("run_ii.mat")

%% Analysis and plot
Analysis_plot
