%% Run the program to get the output
clear
clc

% Path to ''main'' file
oper_filepath = fullfile(fileparts(fileparts(pwd)), 'Main');
% Path back to this file
current_filepath = pwd;
save(fullfile(oper_filepath,'curpath.mat'),'current_filepath')

load('soilprop.mat')
%% scenarios
All_results = cell(1, length(K17(1,1,:)));
save("All_results","All_results")

for ii = 1:length(K17(1,1,:))
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
    save("simulation time.mat",'T_simulation')

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

delete("mydata.mat")
delete("run_ii.mat")

%% Analysis and plot
% Analysis_plot
