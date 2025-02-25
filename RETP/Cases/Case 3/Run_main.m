%% Run the program to get the output
clear
clc

% Path to ''main'' file
oper_filepath = fullfile(fileparts(fileparts(pwd)), 'Main');
% Path back to this file
current_filepath = pwd;
save(fullfile(oper_filepath,'curpath.mat'),'current_filepath')

%% UH, ESH and PH scenarios
Case_suffixes = {'_2DHeter'};
Output_name = {'2DHeter.mat'};
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
    T_simulation = 58000;
    save("simulation time.mat",'T_simulation')

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
    load(Output_name{ii})
    ds = rmfield(ds, {'re_wmass', 'tp_wmass'});
    save(Output_name{ii},'ds')
    
end

delete("Case_name.mat")
delete("run_ii.mat")

%% Analysis and plot
Analysis_plot
