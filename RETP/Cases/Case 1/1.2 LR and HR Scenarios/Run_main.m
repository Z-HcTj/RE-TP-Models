%% Run the program to get the output
clear
clc

% Path to ''main'' file
oper_filepath = fullfile(fileparts(fileparts(fileparts(pwd))), 'Main');
% Path back to this file
current_filepath = pwd;
save(fullfile(oper_filepath,'curpath.mat'),'current_filepath')

%% LR and HR scenarios
Case_suffixes = {'_LR', '_HR'};
Output_name = {'1_LR_results.mat', '2_HR_results.mat'};
save("Case_name","Case_suffixes","Output_name")

for ii = 1:length(Case_suffixes)
    save("run_ii","ii")

    % Input parameters for every scenario
    load('Case_name.mat')
    Casename = ['Inputdata' Case_suffixes{ii} '.m'];
    run(Casename)

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

delete("Case_name.mat")
delete("run_ii.mat")

%% Analysis and plot
Analysis_plot
