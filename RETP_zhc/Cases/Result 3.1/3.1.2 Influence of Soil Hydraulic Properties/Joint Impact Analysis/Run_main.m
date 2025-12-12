%% Run the program to get the output
clear
clc

% Path to ''main'' file
oper_filepath = fullfile(fileparts(fileparts(fileparts(fileparts(pwd)))), 'Main');
% Path back to this file
current_filepath = pwd;
save(fullfile(oper_filepath,'curpath.mat'),'current_filepath')

%% Joint Impact
All_results = cell(1, 1000);
save("All_results","All_results")

for ii = 1:length(All_results)
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
    load('All_results.mat')
    load('run_ii.mat')
    All_results{ii} = load("mydata.mat");
    save("All_results","All_results")

end

delete("run_ii.mat")
delete("mydata.mat")

%% Analysis and plot
Analysis_plot

%% Random Forest Out-of-Bag Permutation Error Analysis
Random_forest