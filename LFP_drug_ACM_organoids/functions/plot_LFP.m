function plot_LFP(agg_file1, agg_file2, wells2process, output_LFP, t_index)
    % plot_LFP: Function to plot Local Field Potentials (LFP) from two aggregated files
    % agg_file1: Path to the first aggregated LFP file
    % agg_file2: Path to the second aggregated LFP file
    % wells2process: Wells to be processed, e.g. [1, 2, 3]
    % output_LFP: Directory to save the output plots
    % t_index: Time index for filtering the data

    % Load the first LFP file
    if exist(agg_file1, 'file')
        load(agg_file1, 'LFP', 't_ds'); % Load LFP data
        LFP1 = LFP;  % Store LFP data from the first file
        t_ds1 = t_ds; % Store time data from the first file
    else
        disp('Error: The first LFP file does not exist.');
        return;
    end
    
    % Load the second LFP file
    if exist(agg_file2, 'file')
        load(agg_file2, 'LFP', 't_ds'); % Load LFP data
        LFP2 = LFP;  % Store LFP data from the second file
        t_ds2 = t_ds; % Store time data from the second file
    else
        disp('Error: The second LFP file does not exist.');
        return;
    end

    % Create output directory if it does not exist
    if ~exist(output_LFP, 'dir')
        mkdir(output_LFP);
    end

    % Process each well
    for i = 1:length(wells2process)
        well_index = wells2process(i);  % Current well index

        % Ensure LFP data exists for the current well
        if well_index <= length(LFP1)
            % Sum LFP data for the first file
            LFP_sum1 = sum(LFP1{well_index}, 2)*10^6; % Convert to µV

            % Ensure LFP data exists for the second file
            if well_index <= length(LFP2)
                % Sum LFP data for the second file
                LFP_sum2 = sum(LFP2{well_index}, 2)*10^6; % Convert to µV
            else
                fprintf('Warning: No data for well %d in the second LFP file.\n', well_index);
                continue;
            end

            % Filter time indices based on the provided time range
            time_index1 = t_ds1 >= 0 & t_ds1 <= t_index;  % Time filter for the first dataset
            time_index2 = t_ds2 >= 0 & t_ds2 <= t_index;  % Time filter for the second dataset

            % Create figure for plotting
            figure('Units', 'Inches', 'Position', [1, 1, 10, 4]);

            % Plot LFP from the first file
            p1 = plot(t_ds1(time_index1), LFP_sum1(time_index1), 'LineWidth', 0.001, 'Color', [24,110,189] / 255); % Line color for LFP1
            hold on;

            % Plot LFP from the second file
            p2 = plot(t_ds2(time_index2), LFP_sum2(time_index2), 'LineWidth', 0.001, 'Color', [210,62,18] / 255); % Line color for LFP2
            
            % Set transparency for the plots
            p1.Color(4) = 0.5; % Set alpha for LFP1
            p2.Color(4) = 0.25; % Set alpha for LFP2
            
            % Label axes
            xlabel('Time (s)', 'FontSize', 12, 'FontName', 'Arial');
            ylabel('Voltage (µV)', 'FontSize', 12, 'FontName', 'Arial');

            % Set axis properties
            ax = gca;
            ax.FontSize = 12; % Axis font size
            set(gca, 'FontName', 'Arial'); % Axis font
            xlim([0 t_index]); % Set x-axis limits
            ylim([-45,30]); % Set y-axis limits
            
            % Save the figure
            saveas(gcf, fullfile(output_LFP, sprintf('well_%d.png', well_index)));

            % Close the figure
            close;
        else
            fprintf('Warning: No data for well %d.\n', well_index);
        end
    end
end
