%% Run the program to get the output
clear
clc

% Path to ''main'' file
oper_filepath = fullfile(fileparts(fileparts(fileparts(pwd))), 'Main');
% Path back to this file
current_filepath = pwd;
save(fullfile(oper_filepath,'curpath.mat'),'current_filepath')

%% Single Impact
Output_name = {'CMF.mat' 'CFM.mat' 'FMC.mat'...
    'Homo.C.mat' 'Homo.M.mat' 'Homo.F.mat'};
key_param = [5.03 4.6753 2.9298 5.03 2.9298 4.6753 2.9298 4.6753 5.03 5.03 5.03 5.03 4.6753 4.6753 4.6753 2.9298 2.9298 2.9298;
    2.896 2.3303 1.4852 2.896 1.4852 2.3303 1.4852 2.3303 2.896 2.896 2.896 2.896 2.3303 2.3303 2.3303 1.4852 1.4852 1.4852
    41.4 35.5 17.4 41.4 17.4 35.5 17.4 35.5 41.1 41.1 41.1 41.1 35.5 35.5 35.5 17.4 17.4 17.4];

save("Case_name","Output_name","key_param")

for ii = 1:length(Output_name)
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
