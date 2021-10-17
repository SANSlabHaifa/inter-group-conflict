function [system1_output, system2_output, trials, fs] = ...
    getSignalsAndEventsHOMER(channel_def_file, nirs_file, signal_type)

% [system1_output, system2_output, trials, fs] = ...
%    getSignalsAndEventsHOMER(channel_def_file, nirs_file, signal_type)
%
% Returns arrays containing oxy-Hb/deoxy-Hb concentrations over time in
% each channel for two "Brite" systems whose outputs were recorded in a
% single file, as well as information regarding the timing of trial beginnings
% and endings.
%
% Inputs:
% =======
% channel_def_file- full path to a .mat file containing a definition of the
% montage used during recording. This file should contain two vectors, 'tx'
% and 'rx', with corresponding elements referring to the source and receiver
% constituting a given channel, respectively
% nirs_file- full path to a .nirs file
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
% trials- an array of structures, in which each structure contains two
% fields; 'onset' and 'end'
% fs- sampling rate, in Hz
%
% Written by Michael Nevat, 2020.

load(channel_def_file)

%%%[pathstr,filename,ext] = fileparts(nirs_file);
%%%clc
%%%disp(['Reading file ', filename, ext, '...'])
disp('Reading data...')

if ~exist('signal_type'), signal_type = 'O2Hb';, end

load(nirs_file, '-mat')
if ~exist('procResult')
    disp('Hb concentration data not found')
    return
end

if strcmp(signal_type, 'O2Hb')
    tmp = procResult.dc(:, 1, :);
else
    tmp = procResult.dc(:, 2, :);
end
if (size(tmp, 3) ~= 2 * length(rx))
    disp('Mismatch between Hb concentration data and channel definitions')
    return
end

system1_output = tmp(:, 1 : length(rx));
system2_output = tmp(:, length(rx) + 1 : end);

if (abs(mean(diff(t) - 0.1)) < abs(mean(diff(t) - 0.02)))
    fs = 10;
else
    fs = 50;
end

trial_onset_cols = [];
trial_end_cols = [];
for c = 1 : length(CondNames)
    if strfind(CondNames{c}, 'Trial Start')
        trial_onset_cols = [trial_onset_cols c];
    elseif strfind(CondNames{c}, 'Trial End')
        trial_end_cols = [trial_end_cols c];
    end
end
tmp = sum(s(:, trial_onset_cols), 2);
trial_onset_samples = find(tmp);
tmp = sum(s(:, trial_end_cols), 2);
trial_end_samples = find(tmp);
if (~isempty(find((trial_end_samples - trial_onset_samples) < 0)) | ...
        ~isempty(find((trial_end_samples(1 : end - 1) - trial_onset_samples(2 : end)) > 0)) | ...
        (length(trial_end_samples) > length(trial_onset_samples)))
    disp('Mismatch between markers indicating beginning of trials and markers indicating end of trials')
    return
end
for t = 1 : length(trial_onset_samples)
    trials(t).onset = trial_onset_samples(t);
end
for t = 1 : length(trial_end_samples)
    trials(t).end = trial_end_samples(t);
end