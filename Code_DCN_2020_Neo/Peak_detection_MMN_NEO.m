%% %%%%% Peak detection and Mismatch negativity%%%%%%%%%% NEO
% step after PCA_channel_selection_NEO and manual selection of fitting
% datasets - PCA_check - used for SbefD and Dev trials

%% Intro
addpath('/server/fo2-13/data/FETAL_USER/juliam/Consciousness/Analysis')

[path_1A]=files_1A_NEO;
[path_1B]=files_1B_NEO;
[path_2A]=files_2A_NEO;
[path_2B]=files_2B_NEO;
        
%type of trigger that whats to be evaluated
trigS='Sbef';
trigDA='DevA';
trigDB='DevB';
trigDC='DevC';
trigDD='DevD';

%time before trigger
tbeft=200; 
%time after trigger
taft=3000; 

set(0,'defaultfigurecolor',[1 1 1])
%% read in results
cd '/server/fo2-13/data/FETAL_USER/juliam/Consciousness/Analysis/Results/NEO'
%read in results from previous analysis
load('Results_1A_SbefD_NEO') %load results
load('Results_1A_DevB_NEO')
load('Results_1B_SbefD_NEO')
load('Results_1B_DevA_NEO')
load('Results_2A_SbefD_NEO')
load('Results_2A_DevD_NEO')
load('Results_2B_SbefD_NEO')
load('Results_2B_DevC_NEO')

sf=610.3516; %sampling frequency, in case not calculated before 
%% evaluate results
%load results
%extract info about top 5 sensor positions
[Top5_1A]=Chan_pos_top5_neo(results_1A_SbefD, results_1A_DevB);
[Top5_1B]=Chan_pos_top5_neo(results_1B_SbefD, results_1B_DevA);
[Top5_2A]=Chan_pos_top5_neo(results_2A_SbefD, results_2A_DevD);
[Top5_2B]=Chan_pos_top5_neo(results_2B_SbefD, results_2B_DevC);
%1. position index for NEO data: All brain activity should be somewhere in
%the middle of the sensor array and not at the very top or bottom, as the
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
%% make selected complonents compatible to format from xls sheet
%1A
PCA_check_1A_SbefD=PCA_check_creation_function(sel_comp_1AS, Neo_pos_index_1AS);
PCA_check_1A_Dev=PCA_check_creation_function(sel_comp_1AD, Neo_pos_index_1AD);
%1B
PCA_check_1B_SbefD=PCA_check_creation_function(sel_comp_1BS, Neo_pos_index_1BS);
PCA_check_1B_Dev=PCA_check_creation_function(sel_comp_1BD, Neo_pos_index_1BD);
%2A
PCA_check_2A_SbefD=PCA_check_creation_function(sel_comp_2AS, Neo_pos_index_2AS);
PCA_check_2A_Dev=PCA_check_creation_function(sel_comp_2AD, Neo_pos_index_2AD);
%2B
PCA_check_2B_SbefD=PCA_check_creation_function(sel_comp_2BS, Neo_pos_index_2BS);
PCA_check_2B_Dev=PCA_check_creation_function(sel_comp_2BD, Neo_pos_index_2BD);
%% check if ther is a suitable component for Std and Dev - if not, set counterpart to 0
[PCA_check_1A_SbefD, PCA_check_1A_Dev]=JM_check_matching_PCA_check(PCA_check_1A_SbefD, PCA_check_1A_Dev);
%1B
[PCA_check_1B_SbefD, PCA_check_1B_Dev]=JM_check_matching_PCA_check(PCA_check_1B_SbefD, PCA_check_1B_Dev);
%2A
[PCA_check_2A_SbefD, PCA_check_2A_Dev]=JM_check_matching_PCA_check(PCA_check_2A_SbefD, PCA_check_2A_Dev);
%2B
[PCA_check_2B_SbefD, PCA_check_2B_Dev]=JM_check_matching_PCA_check(PCA_check_2B_SbefD, PCA_check_2B_Dev);

%%

%%
% %read in data from visual evaluation of PCA results
% PCA_check_1A=readtable('PCA_analysis_check_1A_for_matlab_NEO_adapted.xlsx');
% PCA_check_1B=readtable('PCA_analysis_check_1B_for_matlab_NEO_adapted.xlsx');
% PCA_check_2A=readtable('PCA_analysis_check_2A_for_matlab_NEO_adapted.xlsx');
% PCA_check_2B=readtable('PCA_analysis_check_2B_for_matlab_NEO_adapted.xlsx');
% %%
% %PCA_check specific for SbefD and Dev
% PCA_check_1A_SbefD=table2array(PCA_check_1A(:,2));
% PCA_check_1A_Dev=table2array(PCA_check_1A(:,3));
% PCA_check_1B_SbefD=table2array(PCA_check_1B(:,2));
% PCA_check_1B_Dev=table2array(PCA_check_1B(:,3));
% PCA_check_2A_SbefD=table2array(PCA_check_2A(:,2));
% PCA_check_2A_Dev=table2array(PCA_check_2A(:,3));
% PCA_check_2B_SbefD=table2array(PCA_check_2B(:,2));
% PCA_check_2B_Dev=table2array(PCA_check_2B(:,3));
%% addition: remove noisy datasets
%remove all datasets where artifact blocking threshold is corrected to 2pT
%-> these datasets have a large number of high amplitude noise and will not
%be used for further analysis
nD1A=find(cell2mat(results_1A_SbefD{3,2})==2e-12);
PCA_check_1A_SbefD(nD1A)=0;
PCA_check_1A_Dev(nD1A)=0;

nD1B=find(cell2mat(results_1B_SbefD{3,2})==2e-12);
PCA_check_1B_SbefD(nD1B)=0;
PCA_check_1B_Dev(nD1B)=0;

nD2A=find(cell2mat(results_2A_SbefD{3,2})==2e-12);
PCA_check_2A_SbefD(nD2A)=0;
PCA_check_2A_Dev(nD2A)=0;

nD2B=find(cell2mat(results_2B_SbefD{3,2})==2e-12);
PCA_check_2B_SbefD(nD2B)=0;
PCA_check_2B_Dev(nD2B)=0;

%% save for later use 
cd '/server/fo2-13/data/FETAL_USER/juliam/Consciousness/Analysis/Results/NEO'
%1A
save('PCA_check_1A_SbefD_NEO.mat', 'PCA_check_1A_SbefD')
save('PCA_check_1A_Dev_NEO.mat', 'PCA_check_1A_Dev')
%1B
save('PCA_check_1B_SbefD_NEO.mat', 'PCA_check_1B_SbefD')
save('PCA_check_1B_Dev_NEO.mat', 'PCA_check_1B_Dev')
%2A
save('PCA_check_2A_SbefD_NEO.mat', 'PCA_check_2A_SbefD')
save('PCA_check_2A_Dev_NEO.mat', 'PCA_check_2A_Dev')
%2B
save('PCA_check_2B_SbefD_NEO.mat', 'PCA_check_2B_SbefD')
save('PCA_check_2B_Dev_NEO.mat', 'PCA_check_2B_Dev')
%%  find peak locations - plot 08.08.2018
%1A
JM_peak_location_plot_v2(path_1A, sf, trigS, tbeft, taft, results_1A_SbefD, PCA_check_1A_SbefD)
JM_peak_location_plot_v2(path_1A, sf, trigDB, tbeft, taft, results_1A_DevB, PCA_check_1A_Dev)

%1B
JM_peak_location_plot_v2(path_1B, sf, trigS, tbeft, taft, results_1B_SbefD, PCA_check_1B_SbefD)
JM_peak_location_plot_v2(path_1B, sf, trigDA, tbeft, taft, results_1B_DevA, PCA_check_1B_Dev)

%2A
JM_peak_location_plot_v2(path_2A, sf, trigS, tbeft, taft, results_2A_SbefD, PCA_check_2A_SbefD)
JM_peak_location_plot_v2(path_2A, sf, trigDD, tbeft, taft, results_2A_DevD, PCA_check_2A_Dev)

%2B
JM_peak_location_plot_v2(path_2B, sf, trigS, tbeft, taft, results_2B_SbefD, PCA_check_2B_SbefD)
JM_peak_location_plot_v2(path_2B, sf, trigDC, tbeft, taft, results_2B_DevC, PCA_check_2B_Dev)
%% MM within sequence - tone 3 and 4  - 07.08.2018 - change: 04.09.2018: rms of 5 channels instead of mean - v3
%always deviant - standard
%select result fitting to trigger
%1A
[mm_peaks_1A_SbefD]= MMwithinSequence_v3(path_1A, sf, tbeft, trigS, results_1A_SbefD, PCA_check_1A_SbefD);
[mm_peaks_1A_DevB]= MMwithinSequence_v3(path_1A, sf, tbeft, trigDB, results_1A_DevB, PCA_check_1A_Dev);

%1B
[mm_peaks_1B_SbefD]= MMwithinSequence_v3(path_1B, sf, tbeft, trigS, results_1B_SbefD, PCA_check_1B_SbefD);
[mm_peaks_1B_DevA]= MMwithinSequence_v3(path_1B, sf, tbeft, trigDA, results_1B_DevA, PCA_check_1B_Dev);

%2A
[mm_peaks_2A_SbefD]= MMwithinSequence_v3(path_2A, sf, tbeft, trigS, results_2A_SbefD, PCA_check_2A_SbefD);
[mm_peaks_2A_DevD]= MMwithinSequence_v3(path_2A, sf, tbeft, trigDD, results_2A_DevD, PCA_check_2A_Dev);

%2B
[mm_peaks_2B_SbefD]= MMwithinSequence_v3(path_2B, sf, tbeft, trigS, results_2B_SbefD, PCA_check_2B_SbefD);
[mm_peaks_2B_DevC]= MMwithinSequence_v3(path_2B, sf, tbeft, trigDC, results_2B_DevC, PCA_check_2B_Dev);


cd '/server/fo2-13/data/FETAL_USER/juliam/Consciousness/FAIRY_analysis/Results'
%cd '/server/fo2-13/data/FETAL_USER/juliam/Consciousness/Analysis/Results/NEO'
%save MM peaks
%1A
save('MM_peaks_1A_SbefD_NEO.mat', 'mm_peaks_1A_SbefD')
save('MM_peaks_1A_DevB_NEO.mat', 'mm_peaks_1A_DevB')

%1B
save('MM_peaks_1B_SbefD_NEO.mat', 'mm_peaks_1B_SbefD')
save('MM_peaks_1B_DevA_NEO.mat', 'mm_peaks_1B_DevA')

%2A
save('MM_peaks_2A_SbefD_NEO.mat', 'mm_peaks_2A_SbefD')
save('MM_peaks_2A_DevD_NEO.mat', 'mm_peaks_2A_DevD')

%2B
save('MM_peaks_2B_SbefD_NEO.mat', 'mm_peaks_2B_SbefD')
save('MM_peaks_2B_DevC_NEO.mat', 'mm_peaks_2B_DevC')
%% MM between sequences - 07.08.2018 - change: 04.09.2018: rms of 5 channels instead of mean - v3
%always deviant - standard
%1A
[mm_peaks_1A_between]= MMbetweenSequence_v3(path_1A, sf, tbeft, results_1A_SbefD, results_1A_DevB, PCA_check_1A_SbefD, PCA_check_1A_Dev);

%1B
[mm_peaks_1B_between]= MMbetweenSequence_v3(path_1B, sf, tbeft, results_1B_SbefD, results_1B_DevA, PCA_check_1B_SbefD, PCA_check_1B_Dev);

%2A
[mm_peaks_2A_between]= MMbetweenSequence_v3(path_2A, sf, tbeft, results_2A_SbefD, results_2A_DevD, PCA_check_2A_SbefD, PCA_check_2A_Dev);

%2B
[mm_peaks_2B_between]= MMbetweenSequence_v3(path_2B, sf, tbeft, results_2B_SbefD, results_2B_DevC, PCA_check_2B_SbefD, PCA_check_2B_Dev);

%cd '/server/fo2-13/data/FETAL_USER/juliam/Consciousness/FAIRY_analysis/Results'
cd '/server/fo2-13/data/FETAL_USER/juliam/Consciousness/Analysis/Results/NEO'
%save MM peaks
%1A
save('MM_peaks_1A_between_NEO.mat', 'mm_peaks_1A_between')

%1B
save('MM_peaks_1B_between_NEO.mat', 'mm_peaks_1B_between')

%2A
save('MM_peaks_2A_between_NEO.mat', 'mm_peaks_2A_between')

%2B
save('MM_peaks_2B_between_NEO.mat', 'mm_peaks_2B_between')

%% peak count 24.08.2018
% %1A
% [peak_overview_1A_SbefD]=JM_peak_count_v2(results_1A_SbefD, PCA_check_1A_SbefD);
% [peak_overview_1A_DevB]=JM_peak_count_v2(results_1A_DevB, PCA_check_1A_Dev);
% 
% %1B
% [peak_overview_1B_SbefD]=JM_peak_count_v2(results_1B_SbefD, PCA_check_1B_SbefD);
% [peak_overview_1B_DevA]=JM_peak_count_v2(results_1B_DevA, PCA_check_1B_Dev);
% 
% %2A
% [peak_overview_2A_SbefD]=JM_peak_count_v2(results_2A_SbefD, PCA_check_2A_SbefD);
% [peak_overview_2A_DevD]=JM_peak_count_v2(results_2A_DevD, PCA_check_2A_Dev);
% 
% %2B
% [peak_overview_2B_SbefD]=JM_peak_count_v2(results_2B_SbefD, PCA_check_2B_SbefD);
% [peak_overview_2B_DevC]=JM_peak_count_v2(results_2B_DevC, PCA_check_2A_Dev);
% 
% cd '/server/fo2-13/data/FETAL_USER/juliam/Consciousness/Analysis/Results/NEO'
% %save peak counts
% %1A
% save('Peak_overview_1A_SbefD_NEO.mat', 'peak_overview_1A_SbefD')
% save('Peak_overview_1A_DevB_NEO.mat', 'peak_overview_1A_DevB')
% 
% %1B
% save('Peak_overview_1B_SbefD_NEO.mat', 'peak_overview_1B_SbefD')
% save('Peak_overview_1B_DevA_NEO.mat', 'peak_overview_1B_DevA')
% 
% %2A
% save('Peak_overview_2A_SbefD_NEO.mat', 'peak_overview_2A_SbefD')
% save('Peak_overview_2A_DevD_NEO.mat', 'peak_overview_2A_DevD')
% 
% %2B
% save('Peak_overview_2B_SbefD_NEO.mat', 'peak_overview_2B_SbefD')
% save('Peak_overview_2B_DevC_NEO.mat', 'peak_overview_2B_DevC')

%% 03.09.2018 count all peaks that exist for 3rd and 4th tone to calculate a MM within 
% 
% %1A
% [Counted_MM_peaks_within_1A_SbefD]=JM_peak_count_MM_within(results_1A_SbefD, PCA_check_1A_SbefD);
% [Counted_MM_peaks_within_1A_DevB]=JM_peak_count_MM_within(results_1A_DevB, PCA_check_1A_Dev);
% 
% %1B
% [Counted_MM_peaks_within_1B_SbefD]=JM_peak_count_MM_within(results_1B_SbefD, PCA_check_1B_SbefD);
% [Counted_MM_peaks_within_1B_DevA]=JM_peak_count_MM_within(results_1B_DevA, PCA_check_1B_Dev);
% 
% %2A
% [Counted_MM_peaks_within_2A_SbefD]=JM_peak_count_MM_within(results_2A_SbefD, PCA_check_2A_SbefD);
% [Counted_MM_peaks_within_2A_DevD]=JM_peak_count_MM_within(results_2A_DevD, PCA_check_2A_Dev);
% 
% %2B
% [Counted_MM_peaks_within_2B_SbefD]=JM_peak_count_MM_within(results_2B_SbefD, PCA_check_2B_SbefD);
% [Counted_MM_peaks_within_2B_DevC]=JM_peak_count_MM_within(results_2B_DevC, PCA_check_2B_Dev);
% 
% cd '/server/fo2-13/data/FETAL_USER/juliam/Consciousness/Analysis/Results/NEO'
% %save peak counts
% %1A
% save('Counted_MM_peaks_within_1A_SbefD_NEO.mat', 'Counted_MM_peaks_within_1A_SbefD')
% save('Counted_MM_peaks_within_1A_DevB_NEO.mat', 'Counted_MM_peaks_within_1A_DevB')
% 
% %1B
% save('Counted_MM_peaks_within_1B_SbefD_NEO.mat', 'Counted_MM_peaks_within_1B_SbefD')
% save('Counted_MM_peaks_within_1B_DevA_NEO.mat', 'Counted_MM_peaks_within_1B_DevA')
% 
% %2A
% save('Counted_MM_peaks_within_2A_SbefD_NEO.mat', 'Counted_MM_peaks_within_2A_SbefD')
% save('Counted_MM_peaks_within_2A_DevD_NEO.mat', 'Counted_MM_peaks_within_2A_DevD')
% 
% %2B
% save('Counted_MM_peaks_within_2B_SbefD.mat', 'Counted_MM_peaks_within_2B_SbefD')
% save('Counted_MM_peaks_within_2B_DevC.mat', 'Counted_MM_peaks_within_2B_DevC')
% 
% %% peak count between 03.09.2018
% %1A
% [tone4_SbefD, filename]=JM_peak_count_MM_for_between(results_1A_SbefD, PCA_check_1A_SbefD);
% [tone4_Dev]=JM_peak_count_MM_for_between(results_1A_DevB, PCA_check_1A_Dev);
% 
% Peaks_between=[tone4_SbefD, tone4_Dev];
% %delete any row that is empty for either og the tones
% [idxy, idxx]=find(Peaks_between==0);
% idx=unique(idxy); %only look for y values as whole rows are deleted anyways
% tone4_SbefD(idx,:)=[];
% tone4_Dev(idx,:)=[];
% filename(idx,:)=[];
% %combine findings into table
% Counted_MM_peaks_between_1A=table(filename, tone4_SbefD, tone4_Dev);
% 
% %1B
% [tone4_SbefD, filename]=JM_peak_count_MM_for_between(results_1B_SbefD, PCA_check_1B_SbefD);
% [tone4_Dev]=JM_peak_count_MM_for_between(results_1B_DevA, PCA_check_1B_Dev);
% 
% Peaks_between=[tone4_SbefD, tone4_Dev];
% %delete any row that is empty for either og the tones
% [idxy, idxx]=find(Peaks_between==0);
% idx=unique(idxy); %only look for y values as whole rows are deleted anyways
% tone4_SbefD(idx,:)=[];
% tone4_Dev(idx,:)=[];
% filename(idx,:)=[];
% %combine findings into table
% Counted_MM_peaks_between_1B=table(filename, tone4_SbefD, tone4_Dev);
% 
% %2A
% [tone4_SbefD, filename]=JM_peak_count_MM_for_between(results_2A_SbefD, PCA_check_2A_SbefD);
% [tone4_Dev]=JM_peak_count_MM_for_between(results_2A_DevD, PCA_check_2A_Dev);
% 
% Peaks_between=[tone4_SbefD, tone4_Dev];
% %delete any row that is empty for either og the tones
% [idxy, idxx]=find(Peaks_between==0);
% idx=unique(idxy); %only look for y values as whole rows are deleted anyways
% tone4_SbefD(idx,:)=[];
% tone4_Dev(idx,:)=[];
% filename(idx,:)=[];
% %combine findings into table
% Counted_MM_peaks_between_2A=table(filename, tone4_SbefD, tone4_Dev);
% 
% %2B
% [tone4_SbefD, filename]=JM_peak_count_MM_for_between(results_2B_SbefD, PCA_check_2B_SbefD);
% [tone4_Dev]=JM_peak_count_MM_for_between(results_2B_DevC, PCA_check_2B_Dev);
% 
% Peaks_between=[tone4_SbefD, tone4_Dev];
% %delete any row that is empty for either og the tones
% [idxy, idxx]=find(Peaks_between==0);
% idx=unique(idxy); %only look for y values as whole rows are deleted anyways
% tone4_SbefD(idx,:)=[];
% tone4_Dev(idx,:)=[];
% filename(idx,:)=[];
% %combine findings into table
% Counted_MM_peaks_between_2B=table(filename, tone4_SbefD, tone4_Dev);