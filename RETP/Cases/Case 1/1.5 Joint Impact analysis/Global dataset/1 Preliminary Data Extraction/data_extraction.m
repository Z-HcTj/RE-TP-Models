% This document performs preliminary processing on the global soil 
% hydraulic properties dataset from Gupta et al. (2022), extracting 
% key parameters (including n, α, k, etc.) required for the study.
clear
clc

%% Preliminary Data Extraction
% Read the Excel data
filename = 'WRC_data_gooddata.xlsx'; % Excel file name
sheet = 'Sheet1'; % Worksheet name
opts = detectImportOptions(filename, 'Sheet', sheet);
opts = setvartype(opts, 'layer_id', 'char'); % Specify 'layer_id' as text type
data = readtable(filename, opts);

% Extract the required columns
% Assume the column names are 'layer_id', 'tex_psda', 'alpha', 'n'
columns_to_extract = {'layer_id', 'tex_psda','ksat_lab', 'alpha', 'n'};
extracted_data = data(:, columns_to_extract);

% Remove duplicate rows based on 'layer_id'
[~, unique_indices] = unique(extracted_data.layer_id); % Get indices of unique 'layer_id'
unique_soil_data = extracted_data(unique_indices, :);

% Filtering rule: 'ksat_lab' column should not be 'NA' and should be convertible to a valid number
is_numeric_ksat = ~strcmp(unique_soil_data.ksat_lab, 'NA') & ...
                  ~isnan(str2double(unique_soil_data.ksat_lab));

% Create a new table variable that only contains rows with valid 'ksat_lab' values
unique_soil_data_kinside = unique_soil_data(is_numeric_ksat, :);
unique_soil_data_kinside.ksat_lab = str2double(unique_soil_data_kinside.ksat_lab);

% Save the processed data to a new file for further analysis
output_filename = 'Cleaned_WRC_Data.xlsx';
writetable(unique_soil_data_kinside, output_filename);