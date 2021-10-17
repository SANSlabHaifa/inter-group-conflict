function [system1_output, system2_output, events, fs] = ...
    getSignalsAndEvents(channel_def_file, xls_file, signal_type)

% [system1_output, system2_output, events, fs] = ...
%    getSignalsAndEvents(channel_def_file, xls_file, signal_type)
%
% Returns arrays containing oxy-Hb/deoxy-Hb concentrations over time in
% each channel for two "Brite" systems whose outputs were recorded in a
% single file, as well as information regarding markers that were entered
% during the recording.
%
% Inputs:
% =======
% channel_def_file- full path to a .mat file containing a definition of the
% montage used during recording. This file should contain two vectors, 'tx'
% and 'rx', with corresponding elements referring to the source and receiver
% constituting a given channel, respectively
% xls_file- full path to an MS-Excel file, generated using "Oxysoft" sofware's
% "export" option, and containing measures of oxy-Hb and deoxy-Hb over time
% signal_type (optional, default = 'O2Hb')- a string indicating which
% signals should be extracted. Allowed values- 'O2Hb' (indicating oxy-Hb
% concentrations should be extracted), and 'HHb' (indicating deoxy-Hb
% concentrations should be extracted)
% 
% Outputs:
% ========
% system1_output- array containing values of the requested output for one
% system, with rows corresponding to samples, and columns to channels
% system2_output- array containing values of the requested output for the
% second system, with rows corresponding to samples, and columns to channels
% events- an array of structures (e.g., which has been returned by the
% function 'getSignalsAndEvents'), in which each structure contains two
% fields; 'name', containing the label recorded by the "Oxysoft" software,
% and 'location', containing the index of the corresponding sample
% fs- sampling rate, in Hz
%
% Written by Michael Nevat, 2019

load(channel_def_file)

if ~exist('signal_type'), signal_type = 'O2Hb';, end

[~,Text, all]=  xlsread(xls_file);
[r, c] = find(strcmp(Text, 'Device ids'));
system1 = num2str(cell2mat(all(r, 2)));
system2 = num2str(cell2mat(all(r, 3)));
[r, c] = find(strcmp(Text, 'Data file sample rate'));
fs = cell2mat(all(r, 2));

for r = 1 : size(all, 1)
    temp = [all{r, :}];
    temp(isnan(temp)) = [];
    if temp == 1 : length(temp)
        last_col = length(temp);
        data_row_index = r + 1;
        n_samples = size(all, 1) - r;
        break
    end
end

disp(['Reading ' , signal_type, ' concentrations...'])

for ch = 1 : length(rx)
    disp(['Channel ', num2str(ch)])
    
    str1 = ['[' (system1) '] Rx' num2str(rx(ch)) ' - Tx' num2str(tx(ch)) ' ' signal_type];
    str2 = ['[' (system2) '] Rx' num2str(rx(ch)) ' - Tx' num2str(tx(ch)) ' ' signal_type];
    
    location_str = strfind(Text, str1);
    if sum([location_str{:}]) == 0
        warning (['No data found for system ', system1, ', channel ', num2str(ch)])
        system1_output(:, ch) = zeros(n_samples, 1);
    else
        Index = find(not(cellfun('isempty', location_str)));
        [col, row] = ind2sub(size(Text), Index);
        data_col_index = all{col, row - 1};
        system1_output(:, ch) = ([all{data_row_index : end, data_col_index}]);
    end
    
    location_str = strfind(Text, str2);
    if sum([location_str{:}]) == 0
        warning (['No data found for system ', system2, ', channel ', num2str(ch)])
        system2_output(:, ch) = zeros(n_samples, 1);
    else
        Index = find(not(cellfun('isempty', location_str)));
        [col, row] = ind2sub(size(Text), Index);
        data_col_index = all{col, row - 1};
        system2_output(:, ch) = ([all{data_row_index : end, data_col_index}]);
    end
    
end

events.location = [];
event_col = ([all(data_row_index : end, last_col)]) ;
for c = 1 : length(event_col)
    if (~isnan(event_col{c}) & ~isempty(event_col{c})& isempty(strfind(event_col{c}, 'NULL')))
        events.location = [events.location; c];
    end
end
events.names = event_col(events.location);