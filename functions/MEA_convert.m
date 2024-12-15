function MEA_convert(raw_file, output_folder, wells2process)
%   MEA_convert(raw_file, output_folder, wells2process)
%       raw_file: file path of the .raw file to be unpacked
%       output_folder: file path of the output folder to be created
%       wells2process: subset of wells to unpack
warning('off')

% Set up intake of raw data
if ~exist(output_folder, 'dir')
    % If the folder doesn't exist, create it and do conversion
    mkdir(output_folder);
end

% Load the Axion file into MATLAB
disp('Loading .raw file...')
file = AxisFile(raw_file);

% Compute plate dimensions from file header
num_wells = length(unique([file.DataSets.ChannelArray.Channels.WellRow])) * ...
        length(unique([file.DataSets.ChannelArray.Channels.WellColumn]));    

% Handle unexpected number of wells (commented out)
% if num_wells ~= 12 && num_wells ~= 48
%     % If the number of channels/well combination is strange, default to 12 wells
%     disp('Weird number of channel/well combination')
%     num_wells = 12;
% end

% Parameters & loading data into MATLAB
% if num_wells == 12
%     % Dimensions of the plate (# wells per row, col)
%     plate_dim = [3, 4];
%     % Number of channels per well
%     num_chan = 64;
%     row = {'A', 'B', 'C'};
%     col = {'1', '2', '3', '4'};
% elseif num_wells == 48
%     % Dimensions of the plate (# wells per row, col)
%     plate_dim = [6, 8]; %[6,8] or [3, 4]
%     % Number of channels per well
%     num_chan = 16; % 16 or 64
%     row = {'A', 'B', 'C', 'D', 'E', 'F'};
%     col = {'1', '2', '3', '4', '5', '6', '7', '8'};
% end

if num_wells == 6
    % Dimensions of the plate (# wells per row, col)
    plate_dim = [2, 3];
    % Number of channels per well
    num_chan = 64;
    row = {'A', 'B'};
    col = {'1', '2', '3'};
end

% Gather wells & save into .mat files
for k = 1:plate_dim(1)
    for l = 1:plate_dim(2)
        
        % Extract wells
        disp(fprintf('Well: (%i,%i)', k, l));        
        % Calculate the well number
        well_number = (k - 1) * plate_dim(2) + l; % Calculate index
        outname = fullfile(output_folder, sprintf('well_%d.mat', well_number));
        
        if any(wells2process == ((k - 1) * plate_dim(2) + l))                        
            % If we want to convert this well
            if ~exist(char(outname), 'file')
                
                % Load the well as the dataset is too large
                data = file.DataSets.LoadData([row{k} col{l}]);                                
                if all(all(squeeze([cellfun(@isempty, data(k, l, :, :))])))
                    % If there is nothing in the well, skip to the next iteration
                    disp(["Whole well (", num2str(k), ", ", num2str(l), ") empty, skipped."]);
                    continue
                    
                elseif any(any(squeeze([cellfun(@isempty, data(k, l, :, :))])))
                    % Some channels are empty, have to go through each channel
                    % Fill empty channels with 0s
                    disp('Empty channels encountered. Sorting...')
                    
                    % Find first non-empty channel
                    [r, c] = find(1 - squeeze(cellfun(@isempty, data(k, l, :, :))), 1);
                    
                    % Get time vector and data length from here
                    t = data{k, l, r, c}.GetTimeVector;
                    data_len = length(t);
                    
                    % Pre-initialize MEA array
                    MEA = zeros(data_len, num_chan);
                    chan_dims = size(squeeze(data(k, l, :,:)));                    
                    for rr = 1:chan_dims(1)
                        for cc = 1:chan_dims(2)
                            if ~isempty(data{k, l, rr, cc})
                                % Get data of non-empty channel
                                MEA(:, (cc - 1) * chan_dims(1) + rr) = data{k, l, rr, cc}.GetVoltageVector;
                            else
                                disp(["(", num2str(rr), ",", num2str(cc), ")", " empty."]);
                            end
                        end
                    end                                                   
                else
                    % All channels have data, just convert and save
                    temp = [data{k, l, :, :}];
                    MEA = temp.GetVoltageVector;
                    t = temp(1).GetTimeVector;
                end                
                save(char(outname), 't', 'MEA', '-v7.3')
            else
                disp('Well conversion already exists.');
            end
        else
            disp('Well skipped by user.');
        end                
    end
end
