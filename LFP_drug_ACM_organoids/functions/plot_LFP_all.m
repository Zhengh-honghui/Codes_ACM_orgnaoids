function plot_LFP_all(agg_file, wells2process, output_LFP, raw_file_name)
    % plot_LFP: Function to plot LFP data
    % agg_file: Path to the LFP aggregated file
    % wells2process: List of wells to process, e.g., [1, 2, 3]
    % output_LFP: Directory to save the output plots
    % raw_file_name: Name of the raw file for labeling the plots

    % Specify the file path
    file_path = agg_file; % Path to the aggregated LFP file

    % Check if the file exists
    if exist(file_path, 'file')
        load(file_path, 'LFP', 't_ds'); % Load LFP data and time stamps

        % Create output directory if it does not exist
        if ~exist(output_LFP, 'dir')
            mkdir(output_LFP);
        end
        
        % Create figure for plots
        figure('Units', 'Inches', 'Position', [1, 1, 6, 10]);

        % Loop through each well to plot LFP
        for i = 1:length(wells2process)
            well_index = wells2process(i);  % Current well index
            
            % Check if the well index is within available LFP data
            if well_index <= length(LFP)
                % Sum the LFP data across channels 
                LFP_sum = sum(LFP{well_index}, 2); % Resulting in a 12500x1 vector

                % Create subplot for each well
                subplot(length(wells2process), 1, i); % Create subplot for each well
                
                % Plot LFP data
                plot(t_ds, LFP_sum, 'LineWidth', 0.1, 'Color', [0 0.4470 0.7410]); % Plot LFP
                xlabel('Time (s)', 'FontSize', 12, 'FontName', 'Arial'); % X-axis label
                ylabel('Voltage (µV)', 'FontSize', 12, 'FontName', 'Arial'); % Y-axis label
                
                % Set properties for the axis
                ax = gca; 
                ax.FontSize = 12; % Font size for axis
                set(gca, 'FontName', 'Arial'); % Font for axis
                xlim([min(t_ds) max(t_ds)]); % Set limits for x-axis
                ylim([min(LFP_sum) max(LFP_sum)]); % Set limits for y-axis
                
                % Set title for the current subplot
                title(sprintf('%s_%d', raw_file_name, well_index), 'FontSize', 12, 'FontName', 'Arial');

            else
                fprintf('Warning: Well index %d exceeds available data\n', well_index);
            end
        end

        % Save the figure as a PNG file in the specified output directory
        saveas(gcf, fullfile(output_LFP, sprintf('%s_LFP_All_Wells.png', raw_file_name)));
        
        % Close the figure
        close;
    else
        disp('Error: LFP aggregated file not found');
    end
end
