%% Run the program to get the output
clear
clc

% Path to ''main'' file
oper_filepath = fullfile(fileparts(fileparts(fileparts(fileparts(pwd)))), 'Main');
% Path back to this file
current_filepath = pwd;
save(fullfile(oper_filepath,'curpath.mat'),'current_filepath')

%% Single Impact
Output_name = {'Heter. n=1.5.mat' 'Heter. n=1.7.mat' 'Heter. n=2.mat' 'Heter. n=3.mat'...
    'Heter. n=4.4.mat' 'Heter. n=7.mat'};
key_param = [1.5 1.7 2 3 4.4 7];
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
