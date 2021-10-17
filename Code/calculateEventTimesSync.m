function [blocks, trials] = calculateEventTimesSync(events, log_path, ...
    num_conditions, trials_per_condition, task_start_marker, fs)

% [blocks, trials] = calculateEventTimesSync(events, log_path, num_conditions, ...
%    trials_per_condition, task_start_marker, fs)
%
% Receives an array of structures containing information regarding the markers entered into
% Artinis's "Oxysoft" software during the course of the "Sync" task.
% Reads the log files generated by the task, and determines the indices of samples in the fNIRS
% corresponding the onset and end of each trial and block/condition, which are returned in an 
% array of structures.
%
% Inputs:
% =======
% events- an array of structures (e.g., which has been returned by the
% function 'getSignalsAndEvents'), in which each structure contains two
% fields; 'name', containing the label recorded by the "Oxysoft" software,
% and 'location', containing the index of the corresponding sample.
% log_path- path to the folder containing the log files generated by the
% task.
% num_conditions- number of different blocks/conditions (typically, 3)
% trials_per_condition- number of trials in each condition
% task_start_marker- string containing the character combination
% corresponding to the marker indicating the beginning of the task
% fs (optional)- sampling frequency, in Hz (typically, 10Hz)
%
% Outputs:
% ========
% blocks- an array of structures, in which each structure contains two
% fields; 'onset' and 'end'
% trials- an array of structures, in which each structure contains two
% fields; 'onset' and 'end'
%
% Written by Michael Nevat, 2019.

if (nargin < 5), fs = 10;, end

task_start = ...
    events.location(find(strcmp(events.names, task_start_marker) ...
    + strcmp(events.names, [task_start_marker, ' '])));

for t = 1 : num_conditions * trials_per_condition
    clc
    disp(['Reading log file for trial no. ', num2str(t), '...'])
    f = dir([log_path, '\*_', num2str(t), '.log']);
    fid = fopen([log_path, '\', f(end).name], 'r');
    if (fid < 0)
        disp('Couldn''t find log file')
    else
        if (mod(t, trials_per_condition) == 1)
            for l = 1 : 13
                fgetl(fid);
            end
            tmp = fgetl(fid);
            if (t == 1), task_onset = datevec(tmp(1 : 26));, end
            fgetl(fid);
            tmp = fgetl(fid);
            %%%trials(c, t).onset = round(events.location(1) + etime(datevec(tmp(1 : 26)), task_onset) ...
            trials(t).onset = round(task_start + etime(datevec(tmp(1 : 26)), task_onset) ...
                * fs);
            tmp = fgetl(fid);
            trials(t).end = round(task_start + etime(datevec(tmp(1 : 26)), task_onset) ...
                * fs);
            blocks(ceil(t / trials_per_condition)).onset = trials(t).onset;
        else
            for l = 1 : 14
                fgetl(fid);
            end
            tmp = fgetl(fid);
            trials(t).onset = round(task_start + etime(datevec(tmp(1 : 26)), task_onset) ...
                * fs);
            tmp = fgetl(fid);
            trials(t).end = round(task_start + etime(datevec(tmp(1 : 26)), task_onset) ...
                * fs);
            if (mod(t, trials_per_condition) == 0), blocks(ceil(t / trials_per_condition)).end = trials(t).end;, end
        end
        fclose(fid);
    end
end