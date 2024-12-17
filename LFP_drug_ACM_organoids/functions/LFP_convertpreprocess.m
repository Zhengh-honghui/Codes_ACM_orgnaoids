function LFP_convertpreprocess(batch_folder, output_LFP)
    % LFP_convertpreprocess: Preprocess and convert RAW data files
    % batch_folder: Path to the folder containing .raw files
    % wells2process: Wells to be processed
    
    % List .raw files in the batch directory
    raw_files_first5 = dir(fullfile(batch_folder, '1209first5_(*)(*).raw'));
    raw_files_last5 = dir(fullfile(batch_folder, '1209last5_(*)(*).raw'));

    % Combine the lists of raw files
    raw_files = [raw_files_first5; raw_files_last5];

    wells2process = [1, 2, 3, 4, 5, 6];  % Specify wells to process
    % Check if any .raw files were found
    if isempty(raw_files)
        disp('Error: No .raw files found.');
        return;
    end
    
    % Process each .raw file
    for i = 1:length(raw_files)
        raw_file = fullfile(batch_folder, raw_files(i).name); % Full path to the current .raw file
        
        % Create an output folder for the raw file
        output_folder = fullfile(batch_folder, strtok(raw_files(i).name, '.'));
        
        % Display the current processing file
        disp(['Processing: ', raw_file]);
        
        % Create output directory if it does not exist
        if ~exist(output_folder, 'dir')
            mkdir(output_folder);
        end
        
        % Convert the raw file using MEA_convert
        MEA_convert(raw_file, output_folder, wells2process);
        data_folder = output_folder;
        
        % Process the converted data using MEA_process
        MEA_process(data_folder, wells2process);
        
        % Specify the location of the aggregated file
        agg_file = fullfile(data_folder, 'LFP_Sp.mat');
        
        % Plot the LFP data
        % Save plots in the output folder under LFP_Well_X
        plot_LFP_all(agg_file, wells2process, output_LFP, strtok(raw_files(i).name, '.'));
    end
    
    disp('Processing completed.');
end
