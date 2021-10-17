function [Cunweighted, Cweighted] = ...
    calculateMeanCoherencePerTrial(channel_def_file, X1_O, X2_O, ...
    trials1, trials2, channels, samples1, samples2, min_period, max_period, ...
    disp_str, fs1, fs2)

% [Cunweighted, Cweighted] = ...
%     calculateMeanCoherencePerTrial(channel_def_file, X1_O, X2_O, ...
%     trials1, trials2, channels, samples1, samples2, min_period, max_period, ...
%     disp_str, fs1, fs2)
%
% Calculates the mean coherence throughout the duration of each trial, in a
% specified frequency band, for each combination of channels, disregarding
% combinations containing noisy channels.
%
% Inputs:
% =======
% channel_def_file- full path to a .mat file containing a definition of the
% montage used during recording. This file should contain two vectors, 'tx'
% and 'rx', with corresponding elements referring to the source and receiver
% constituting a given channel, respectively
% X1_O- a two-dimensional array containing oxy-Hb concentrations measured
% by the 1st fNIRS system. Rows represent samples, and columns represent channels
% X2_O- a two-dimensional array containing oxy-Hb concentrations measured
% by the 2nd fNIRS system. Rows represent samples, and columns represent channels
% trials1- an array of structures containing the indices of fNIRS samples
% corresponding the onset and end of each trial in the 1st system
% trials2- an array of structures containing the indices of fNIRS samples
% corresponding the onset and end of each trial in the 2nd system
% channels- binary vector indicating which channels should be included in the analyses
% samples1- binary vector indicating which samples from the 1st system
% should be included in the analyses
% samples2- binary vector indicating which samples from the 2nd system
% should be included in the analyses
% min_period- minimal period, in seconds, for calculation of average
% coherence coefficients
% max_period- maximal period, in seconds, for calculation of average
% coherence coefficients
% disp_str- a string to be displayed on the screen during each iteration 
% fs1 (optional)- sampling rate for 1st system, in Hz (typically, 10Hz)
% fs2 (optional)- sampling rate for 2nd system, in Hz (typically, 10Hz)
%
% Outputs:
% ========
% Cunweighted- a vector containing mean coherence values for each trial in each 
% channel combination (i.e., combination 1 trial 1, combination 1 trial 2, ...
% combination n trial 1, combination n trial 2, ...). Mean values are
% calculated across valid coefficients in period within the specified
% range, and then the (unweighted) average is calculated across periods
% Cweighted- a vector containing mean coherence values for each trial in each 
% channel combination (i.e., combination 1 trial 1, combination 1 trial 2, ...
% combination n trial 1, combination n trial 2, ...). Mean values are
% calculated across valid coefficients in period within the specified
% range, and then the (weighted) average is calculated across periods
%
% Written by Michael Nevat, 2021

load(channel_def_file)

%%% if (length(trials1) ~= length(trials2))
%%%    disp('Mismatch in number of trials. Aborting')
%%%    return
%%% end

addpath('.\wtc-r16')
%%%addpath('..\Preprocessing')
if (nargin < 12), fs2 = 10;, end
if (nargin < 11), fs1 = 10;, end

if (size(samples1, 1) == 1), samples1 = samples1'; end
if (size(samples2, 1) == 1), samples2 = samples2'; end
if (size(samples1 , 1) < size(X1_O, 1))
    short_vector = samples1;
    samples1 = ones(size(X1_O, 1), 1);
    samples1(1 : length(short_vector)) = short_vector;
elseif (size(samples1 , 1) > size(X1_O, 1))
    samples1 = samples1(1 : size(X1_O, 1));
end
if (size(samples2 , 1) < size(X2_O, 1))
    short_vector = samples2;
    samples2 = ones(size(X2_O, 1), 1);
    samples2(1 : length(short_vector)) = short_vector;
elseif (size(samples2 , 1) > size(X2_O, 1))
    samples2 = samples2(1 : size(X2_O, 1));
end
S1 = [X1_O samples1];
S2 = [X2_O samples2];
[S1, S2, trials, fs] = alignInputs(S1, S2, trials1, trials2, fs1, fs2);
X1_O = S1(:, 1 : end - 1);
X2_O = S2(:, 1 : end - 1);
samples = mean([S1(:, end) S2(:, end)], 2)';
samples(find(samples < 1)) = NaN;

num_trials = length(trials);
if (num_trials == 1)
    tmp = cell2mat(struct2cell(trials));
    mean_trial_duration = (tmp(2) - tmp(1)) / fs;
else
    tmp = cell2mat(permute(struct2cell(trials), [2 3 1]));
    trial_durations = (tmp(:, :, 2) - tmp(:, :, 1)) / fs;
    mean_trial_duration = mean(mean(trial_durations));
end

Cunweighted = 9999 * ones(1, length(rx) ^ 2 * num_trials);
Cweighted = Cunweighted;

if (length(channels) ~= 2 * length(tx))
    disp('Problem encountered in specification of noisy channels. Aborting')
    return
end

if (size(channels, 2) == 1), channels = channels'; end
bad_channels = [1 - channels(1 : length(tx)); 1 - channels(length(tx) + 1 : end)];

for c1 = 1 : length(rx)
    x1 = X1_O(:, c1);
    for c2 = 1 : length(rx)
        clc
        disp(disp_str)
        disp('Calculating coherence coefficients...')
        disp(['Channel combination ', num2str((c1 - 1) * length(rx) + c2), '...'])
        x2 = X2_O(:, c2);
        
        t = 0 : 1 / fs : (length(x1) - 1) / fs;
        
        if (~bad_channels(1, c1) & ~bad_channels(2, c2) & ...
                ~isempty(find(X1_O(:, c1))) & ~isempty(find(X2_O(:, c2))))
            try
                [Rsq,period,scale,coi,sig95]=my_wtc([t; x1'], [t; x2'], 'mcc', 0, 'MakeFigure', 0, ...
                    'Dj', 1 / 12, 'S0', 0.2);
                min_sc_tr = min(find(period >= min_period));
                max_sc_tr = max(find(period <= max_period));
                % Discard samples within the COI of samples that have been
                % marked for exclusion
                influenced_samples = round(period * fs) - 1;
                l_coi = ones(length(influenced_samples), max(influenced_samples));
                u_coi = l_coi;
                for p = 1 : size(l_coi, 1)
                    l_coi(p, 1 : influenced_samples(p)) = NaN;
                    u_coi(p, end - influenced_samples(p) + 1 : end) = NaN;
                end
                valid_samples_start = find(diff([1 isnan(samples)]) == -1);
                valid_samples_end = find(diff([isnan(samples) 1]) == 1);
                SM = repmat(samples, length(period), 1);
                for n = 1 : length(valid_samples_start)
                    sample_range = valid_samples_start(n) : ...
                        min(size(SM, 2), valid_samples_start(n) + size(l_coi, 2) - 1);
                    SM(:, sample_range) = ...
                        SM(:, sample_range) .* l_coi(:, 1 : length(sample_range));
                end
                for n = 1 : length(valid_samples_end)
                    sample_range = ...
                        max(1, valid_samples_end(n) - size(l_coi, 2) + 1) : ...
                        valid_samples_end(n);
                    SM(:, sample_range) = ...
                        SM(:, sample_range) .* ...
                        u_coi(:, end - length(sample_range) + 1 : end);
                end
                Rsq = Rsq .* SM;
                                
                for tr = 1 : num_trials
                    % Extract coherence coefficients for the trial, within
                    % the specified band
                    Rsq_trial = ...
                        Rsq(min_sc_tr : max_sc_tr, trials(tr).onset : trials(tr).end);
                    ns = min(size(Rsq_trial, 2), size(l_coi, 2));
                    % Optional- discard coefficients that are influenced by
                    % samples preceding the beginning of the trial
                    Rsq_trial(:, 1 : ns) = ...
                        Rsq_trial(:, 1 : ns) .* ...
                        l_coi(min_sc_tr : max_sc_tr, 1 : ns);
                    
                 
                    % Calculate unweighted mean
                    Cunweighted(((c1 - 1) * length(rx) + c2 - 1) * num_trials + tr) = ...
                        nanmean(nanmean(Rsq_trial, 2));
                    % Calculate weighted mean
                    for sc = 1 : size(Rsq_trial, 1)
                        c_count(sc) = length(find(~isnan(Rsq_trial(sc, :))));
                        c_sum(sc) = ...
                            sum(Rsq_trial(sc, find(~isnan(Rsq_trial(sc, :)))));
                    end
                    Cweighted(((c1 - 1) * length(rx) + c2 - 1) * num_trials + tr) = ...
                        nanmean(c_sum ./ sum(c_count));
                end
            end
        end
        
    end
end