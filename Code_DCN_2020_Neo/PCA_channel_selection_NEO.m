%% 06.11.2018 - NEO -> update 11.02.2019
% only hilbert preprocessing - no comparison with FAIRY at the moment
% processing selectively for each condition (1A, 1B, 2A, 2B)
%list with files that loads all necessary filenames and paths

addpath('/server/fo2-13/data/FETAL_USER/juliam/Consciousness/Analysis')
%use functions from work for EMBC paper
addpath('/server/fo2-13/data/FETAL_USER/juliam/fMEG/EMBC_Analysis')
addpath('/server/fo2-13/data/FETAL_USER/juliam/fMEG/EMBC_Analysis/FirstAttempt')


%files to process
[path_1A, filename_H_1A]=files_1A_NEO;
[path_1B, filename_H_1B]=files_1B_NEO;
[path_2A, filename_H_2A]=files_2A_NEO;
[path_2B, filename_H_2B]=files_2B_NEO;
        
%type of trigger that whats to be evaluated
trigS='Sbef';
trigDA='DevA';
trigDB='DevB';
trigDC='DevC';
trigDD='DevD';

%time before trigger
tbeft=200; 
%time after trigger
taft=3000; %200ms longer then fetal analysis

set(0,'defaultfigurecolor',[1 1 1])
%% PCA channel selection
%1A
[results_1A_SbefD]=JM_PCA_whole_recording_and_SbefD_NEO(path_1A, filename_H_1A, trigS, tbeft, taft);
[results_1A_DevB]=JM_PCA_deviant_trials_NEO(path_1A, filename_H_1A, trigDB, tbeft, taft);

cd '/server/fo2-13/data/FETAL_USER/juliam/Consciousness/Analysis/Results/NEO'
%save results

%1A
save('Results_1A_SbefD_NEO.mat', 'results_1A_SbefD')
save('Results_1A_DevB_NEO.mat', 'results_1A_DevB')

%1B
[results_1B_SbefD]=JM_PCA_whole_recording_and_SbefD_NEO(path_1B, filename_H_1B, trigS, tbeft, taft)
[results_1B_DevA]=JM_PCA_deviant_trials_NEO(path_1B, filename_H_1B, trigDA, tbeft, taft);

cd '/server/fo2-13/data/FETAL_USER/juliam/Consciousness/Analysis/Results/NEO'
%save results

%1B
save('Results_1B_SbefD_NEO.mat', 'results_1B_SbefD')
save('Results_1B_DevA_NEO.mat', 'results_1B_DevA')

%2A
[results_2A_SbefD]=JM_PCA_whole_recording_and_SbefD_NEO(path_2A, filename_H_2A, trigS, tbeft, taft)
[results_2A_DevD]=JM_PCA_deviant_trials_NEO(path_2A, filename_H_2A, trigDD, tbeft, taft);

cd '/server/fo2-13/data/FETAL_USER/juliam/Consciousness/Analysis/Results/NEO'
%save results

%2A
save('Results_2A_SbefD_NEO.mat', 'results_2A_SbefD')
save('Results_2A_DevD_NEO.mat', 'results_2A_DevD')

%2B
[results_2B_SbefD]=JM_PCA_whole_recording_and_SbefD_NEO(path_2B, filename_H_2B,  trigS, tbeft, taft)
[results_2B_DevC]=JM_PCA_deviant_trials_NEO(path_2B, filename_H_2B, trigDC, tbeft, taft);


cd '/server/fo2-13/data/FETAL_USER/juliam/Consciousness/Analysis/Results/NEO'
%save results

%2B
save('Results_2B_SbefD_NEO.mat', 'results_2B_SbefD')
save('Results_2B_DevC_NEO.mat', 'results_2B_DevC')
%% evaluate results
%load results
%extract info about top 5 sensor positions
[Top5_1A]=Chan_pos_top5_neo(Fresults_1A_SbefD, Fresults_1A_DevB);
[Top5_1B]=Chan_pos_top5_neo(Fresults_1B_SbefD, Fresults_1B_DevA);
[Top5_2A]=Chan_pos_top5_neo(Fresults_2A_SbefD, Fresults_2A_DevD);
[Top5_2B]=Chan_pos_top5_neo(Fresults_2B_SbefD, Fresults_2B_DevC);
%1. position index for NEO data: All brain activity should be somewhere in
%the middle of the senso array and not at the very top or bottom, as the
%child is positioned in the craddle accordingly
%"good" components ar elabled with a 1
[Neo_pos_index_1AS, Neo_pos_index_1AD]=Index_neo_position(Top5_1A);
[Neo_pos_index_1BS, Neo_pos_index_1BD]=Index_neo_position(Top5_1B);
[Neo_pos_index_2AS, Neo_pos_index_2AD]=Index_neo_position(Top5_2A);
[Neo_pos_index_2BS, Neo_pos_index_2BD]=Index_neo_position(Top5_2B);
%2. distance within channel clusters
[within_dist_index_1AS, within_dist_index_1AD]=Index_distance_within_chan_clusters_neo(Top5_1A);
[within_dist_index_1BS, within_dist_index_1BD]=Index_distance_within_chan_clusters_neo(Top5_1B);
[within_dist_index_2AS, within_dist_index_2AD]=Index_distance_within_chan_clusters_neo(Top5_2A);
[within_dist_index_2BS, within_dist_index_2BD]=Index_distance_within_chan_clusters_neo(Top5_2B);
%3. overlap with heart activity

% [index_overlap_heart_1AS, index_overlap_heart_1AD]=Index_overlap_heart_neo(Top5_1A);
% [index_overlap_heart_1BS, index_overlap_heart_1BD]=Index_overlap_heart_neo(Top5_1B);
% [index_overlap_heart_2AS, index_overlap_heart_2AD]=Index_overlap_heart_neo(Top5_2A);
% [index_overlap_heart_2BS, index_overlap_heart_2BD]=Index_overlap_heart_neo(Top5_2B);
%4. distance of clusters for standard and deviant
%% selection of 5 sensors taken for evaluation
%% selection 
%1A standard
sel=Neo_pos_index_1AS+within_dist_index_1AS;
[x,y]=find(sel==2);
[~,yy]=unique(y); %if there is more than  selection, take the first one as it has a higher explained variance
sel_comp_1AS=[y(yy)';x(yy)']; %top row: dataset index, bottom row: component (1-3)
%1A deviant
sel=Neo_pos_index_1AD+within_dist_index_1AD;
[x,y]=find(sel==2);
[~,yy]=unique(y); %if there is more than  selection, take the first one as it has a higher explained variance
sel_comp_1AD=[y(yy)';x(yy)']; 

%1B standard
sel=Neo_pos_index_1BS+within_dist_index_1BS;
[x,y]=find(sel==2);
[~,yy]=unique(y); %if there is more than  selection, take the first one as it has a higher explained variance
sel_comp_1BS=[y(yy)';x(yy)']; %top row: dataset index, bottom row: component (1-3)
%1B deviant
sel=Neo_pos_index_1BD+within_dist_index_1BD;
[x,y]=find(sel==2);
[~,yy]=unique(y); %if there is more than  selection, take the first one as it has a higher explained variance
sel_comp_1BD=[y(yy)';x(yy)']; 

%2A standard
sel=Neo_pos_index_2AS+within_dist_index_2AS;
[x,y]=find(sel==2);
[~,yy]=unique(y); %if there is more than  selection, take the first one as it has a higher explained variance
sel_comp_2AS=[y(yy)';x(yy)']; %top row: dataset index, bottom row: component (1-3)
%2A deviant
sel=Neo_pos_index_2AD+within_dist_index_2AD;
[x,y]=find(sel==2);
[~,yy]=unique(y); %if there is more than  selection, take the first one as it has a higher explained variance
sel_comp_2AD=[y(yy)';x(yy)']; 

%2B standard
sel=Neo_pos_index_2BS+within_dist_index_2BS;
[x,y]=find(sel==2);
[~,yy]=unique(y); %if there is more than  selection, take the first one as it has a higher explained variance
sel_comp_2BS=[y(yy)';x(yy)']; %top row: dataset index, bottom row: component (1-3)
%2B deviant
sel=Neo_pos_index_2BD+within_dist_index_2BD;
[x,y]=find(sel==2);
[~,yy]=unique(y); %if there is more than  selection, take the first one as it has a higher explained variance
sel_comp_2BD=[y(yy)';x(yy)']; 