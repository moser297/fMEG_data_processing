

function [results]=JM_PCA_deviant_trials_NEO(path, filename, trig, tbeft, taft)
% function to perform PCA on deviant rtrials as analysis on whole recording
% is only needed once and is done together with standard before deviant
% trials

for h=1:size(path,1)
 currentpath=cell2mat(path(h,:));
 if iscell(filename)
    currentfile=cell2mat(filename(h,:));
 else
     currentfile=filename;
 end
 
 %read in files
[channelinfo, data_new_filtered, ~,sf]=JM_hilbert_read_NEO_E_AB(currentpath, currentfile); % Hilbert
%peak search on trial data
[whole_trial_m, remaining_trials, high5_m, high5_m2, high5_m3, explained, ...
    high5, high52, high53, pks, ptm, pks2, ptm2, pks3, ptm3, ...
    pks_nn, ptm_nn, pks2_nn, ptm2_nn, pks3_nn, ptm3_nn]=JM_PCA_peak_search_triggeravg_NEO_AB(data_new_filtered,...
    channelinfo, sf, currentpath, currentfile, trig, tbeft, taft);
channel_info{h}=channelinfo;
r_trials{h}=remaining_trials;
explained_var{h}=explained;
trial_data_all{h}=whole_trial_m;
trial_data_high5{h}=high5_m;
trial_data_high52{h}=high5_m2;
trial_data_high53{h}=high5_m3;
sensors_high5{h}=high5;
sensors_high52{h}=high52;
sensors_high53{h}=high53;
peaks{h}=pks;
peak_time{h}=ptm;
peaks2{h}=pks2;
peak_time2{h}=ptm2;
peaks3{h}=pks3;
peak_time3{h}=ptm3;
peaks_plusminus{h}=pks_nn;
peak_time_plusminus{h}=ptm_nn;
peaks_plusminus2{h}=pks2_nn;
peak_time_plusminus2{h}=ptm2_nn;
peaks_plusminus3{h}=pks3_nn;
peak_time_plusminus3{h}=ptm3_nn;
end


% search peaks to include into results - peaks later used for plot
[selected_peaks]=JM_FindPeaks(sf, peak_time, peaks);  
[selected_peaks2]=JM_FindPeaks(sf, peak_time2, peaks2); 
[selected_peaks3]=JM_FindPeaks(sf, peak_time3, peaks3); 

[selected_peaks_plusminus]=JM_FindPeaks(sf, peak_time_plusminus, peaks_plusminus);  
[selected_peaks_plusminus2]=JM_FindPeaks(sf, peak_time_plusminus2, peaks_plusminus2); 
[selected_peaks_plusminus3]=JM_FindPeaks(sf, peak_time_plusminus3, peaks_plusminus3);

results={'filename', filename; 'channel_info', channel_info; 'remaining_trials', r_trials'; ...
    'explained_variance', explained_var';  'trial_data', trial_data_all';...
    'rms_PC1', trial_data_high5'; 'rms_PC2', trial_data_high52'; 'rms_PC3', trial_data_high53'; 'sensors_PC1', sensors_high5';...
    'sensors_PC2', sensors_high52'; 'sensors_PC3', sensors_high53';'peak_amplitude_PC1', peaks'; 'peak_time_PC1', peak_time';...
    'peak_amplitude_PC2', peaks2'; 'peak_time_PC2', peak_time2';  'peak_amplitude_PC3', peaks3'; 'peak_time_PC3', peak_time3';...
    'peak_amplitude_plusminus_PC1', peaks_plusminus'; 'peak_time_plusminus_PC1', peak_time_plusminus';...
    'peak_amplitude_plusminus_PC2', peaks_plusminus2'; 'peak_time_plusminus_PC2', peak_time_plusminus2';...
    'peak_amplitude_plusminus_PC3', peaks_plusminus3'; 'peak_time_plusminus_PC3', peak_time_plusminus3';...
    'selected_peaks_PC1', selected_peaks'; 'selected_peaks_PC2', selected_peaks2'; 'selected_peaks_PC3', selected_peaks3';...
    'selected_peaks_plusminus_PC1', selected_peaks_plusminus'; 'selected_peaks_plusminus_PC2', selected_peaks_plusminus2';...
    'selected_peaks_plusminus_PC3', selected_peaks_plusminus3'};


end

 