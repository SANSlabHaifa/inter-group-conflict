function S = readData2Struct...
    (channel_def_file, log_file_base_path, fnirs_data_path, ...
    channel_exclusion_data_path, system_assignment_def_file, file_list, ...
    num_conditions, trials_per_condition, task_start_marker)

% S = readData2Struct...
%    (channel_def_file, log_file_base_path, fnirs_data_path, ...
%    channel_exclusion_data_path, system_assignment_def_file, file_list, ...
%    num_conditions, trials_per_condition, task_start_marker)
%
% Reads fNIRS data (O2Hb concentrations, HHb concentrations, and sampling
% rates), as well as task-related information (onsets and endings of
% trials), and enters the information into an array of structures consisting
% of separate structures for each participant
%
% Inputs:
% =======
% channel_def_file- full path to a .mat file containing a definition of the
% montage used during recording. This file should contain two vectors, 'tx'
% and 'rx', with corresponding elements referring to the source and receiver
% constituting a given channel, respectively
% log_file_base_path- path to a folder containing the folders with the "log" 
% files generated by the task, specified in the "file list" (see below)
% fnirs_data_path- path to a folder containing files with fNIRS data
% (MS-Excel files generated by Artinis's "OxySoft" software and, optionally,
% .nirs files)
% channel_exclusion_data_path- path to a folder containing the MS-Excel files
% specified in the 4th column of the "file list" (see below), which contain
% information regarding channels that are to be excluded from analyses. If
% such files do not exist, pass an empty string (i.e., '') to the function
% system_assignment_def_file- full path to an MS-Excel or .csv file
% specifying the fNIRS system assigned to each participant. The should
% consist of three columns, with the first containing a numerical ID of
% group/dyad, the second containing a numerical identifier of the
% individual participant, and the third containing the serial number of the
% fNIRS, as it appears in the MS-Excel files generated by Artinis's
% "OxySoft" software. The columns should be labelled "Group", "ID", and
% "System", respectively
% file_list- a .csv file, structured as follows:
% * First column- a list of MS-Excel files generated by Artinis' "OxySoft"
% software
% * Second column- a list of folders containing the "log" files generated
% by the task for each group/dyad. The names of these folders should
% contain (among other things) a numerical string identifying the
% group/dyad, which corresponds to records in the file specifying the
% assignment of fNIRS systems (see above)
% * Third column (optional)- a list of .nirs files.
% * Fourth column (optional)- a list of MS-Excel or .csv files containing
% information regarding the channels to be excluded from analyses. These
% files should contain a binary vector, with "1" indicated that a given
% channel is to be excluded from analyses, and "0" indicating it is to be
% included in analyses. The order of the channels should correspond to the
% order of the labels assigned to the channels in the HOMER2 toolbox's GUI
% (e.g., "A1", "B1", etc.). Labels can be included in the files in order
% to ensure correct identification
% num_conditions- number of different blocks/conditions (typically, 3)
% trials_per_condition- number of trials in each condition
% task_start_marker (optional)- string containing the character combination
% corresponding to the marker indicating the beginning of the task.
% This argument should be omitted or left empty (i.e., '') if fNIRS data
% are to be read from the .nirs files specified in the 3rd column of the 
% "file list" (see above)
%
% Outputs:
% S- an array of structures containing fNIRS and task-related information
%
% Written by Michael Nevat, 2021

if (nargin < 9), task_start_marker = '';, end
if isempty(task_start_marker)
    data_source_nirs = 1;
else
    data_source_nirs = 0;
end    

[~, files_n_folders,~] = xlsread(file_list);
[values, labels, raw] = xlsread(system_assignment_def_file);
xls_files = files_n_folders(:, 1);
log_file_folders = files_n_folders(:, 2);
if (size(files_n_folders, 2) >= 3)
    nirs_files = files_n_folders(:, 3);
else
    if data_source_nirs
        disp('Couldn''t read fNIRS data from .nirs files, because no files were specified. Aborting')
        return
    end
    nirs_files = '';
end
if (size(files_n_folders, 2) == 4)
    channel_exclusion_files = files_n_folders(:, 4);
else
    channel_exclusion_files = '';
end

groups = values(:, 1);
participant_ids = values(:, 2);
system_ids = values(:, 3);

S = struct();
p = 0;

for g = 1 : size(files_n_folders, 1)
    
    
    tmp = log_file_folders{g};
    current_group = str2num(tmp(regexp(tmp, '\d')));
    disp(['Group ', num2str(current_group)])
    
    xls_file = [fnirs_data_path, filesep, xls_files{g}];
    [~, ~, raw] = xlsread(xls_file, 'A1:C19');
    [r, ~] = find(strcmp(raw, 'Device ids'));
    systems = [];
    systems(1) = raw{r, 2};
    if ~isnan(raw{r, 3}), systems(2) = raw{r, 3};, end
    
    participants = [];
    participants(1) = ...
        participant_ids(find((groups == current_group) & ...
        (system_ids == systems(1))));
    if (length(systems) == 2)
        participants(2) = ...
            participant_ids(find((groups == current_group) & ...
            (system_ids == systems(2))));
    
    end
    
    o2hb = [];
    hhb = [];
    if data_source_nirs
        nirs_file = [fnirs_data_path, filesep, nirs_files{g}];
        [o2hb(:, :, 1), o2hb(:, :, 2), trials, fs] = ...
            getSignalsAndEventsHOMER(channel_def_file, nirs_file, 'O2Hb');
        [hhb(:, :, 1), hhb(:, :, 2), ~, ~] = ...
            getSignalsAndEventsHOMER(channel_def_file, nirs_file, 'HHb');
    else
        [o2hb(:, :, 1), o2hb(:, :, 2), events, fs] = ...
            getSignalsAndEvents(channel_def_file, xls_file, 'O2Hb');
        [hhb(:, :, 1), hhb(:, :, 2), ~, ~] = ...
            getSignalsAndEvents(channel_def_file, xls_file, 'HHb');
        [blocks, trials] = ...
            calculateEventTimesSync(events, ...
            [log_file_base_path, filesep, log_file_folders{g}], num_conditions, ...
            trials_per_condition, task_start_marker, fs);
    end
    
    clear tIncMan
    if ~isempty(nirs_files)
        load([fnirs_data_path, filesep, nirs_files{g}], '-mat')
    end
    if ~exist('tIncMan')
        tIncMan = ones(size(o2hb, 1), 1);
    elseif (length(tIncMan) < size(o2hb, 1))
        old_length = length(tIncMan);
        new_length = size(o2hb, 1);
        tIncMan(old_length + 1 : new_length) = 1;
    elseif (length(tIncMan) > size(o2hb, 1))
        tIncMan = tIncMan(1 : size(o2hb, 1));
    end
    
    n_channels = size(o2hb, 2);
    if ~isempty(channel_exclusion_files)
        [excluded_channels, ~, ~] = ...
            xlsread([channel_exclusion_data_path, filesep, channel_exclusion_files{g}]);
    else
        excluded_channels = zeros(2 * n_channels, 1);
    end
        
    for s = 1 : length(systems)
        p = p + 1;
        S(p).Group = current_group;
        S(p).Participant = participants(s);
        S(p).SystemID = systems(s);
        S(p).fs = fs;
        S(p).O2Hb = o2hb(:, :, s);
        S(p).HHb = hhb(:, :, s);
        S(p).Conditions = num_conditions;
        S(p).TrialsPerCondition = trials_per_condition;
        S(p).Trials = trials;
        S(p).Samples = tIncMan;
        S(p).Channels = ...
            1 - excluded_channels((s - 1) * n_channels + 1 : s * n_channels);
    end
    
end