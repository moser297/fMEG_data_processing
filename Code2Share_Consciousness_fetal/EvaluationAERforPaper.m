analysis_dir='/server/fo2-13/data/FETAL_USER/juliam/Consciousness/Analysis/NewFAIRYanalysis';
cd(analysis_dir)
%color definitions
set(0,'defaultfigurecolor',[1 1 1])
%standard
s=[0.1 0.1 0.8];
%local & globel
lg=[0.9 0 0.1];
%local
l=[0 0.7 0.8];
%global
g=[1 0.6 0.1];
%local rule
lr=[0 0.8 0];
%local standard
ls=[0.65 0.7 0.6];
%global rule
gr=[0.6 0 0.6];
%global standard
gs=[0.6 0.5 0.6];
grey=[0.4, 0.4, 0.4];
addpath('/server/fo2-13/data/FETAL_USER/juliam/fMEG/AnalysisTools/altmany-export_fig-2763b78')
addpath('/server/fo2-13/data/FETAL_USER/juliam/Consciousness/Analysis')

%%

sf=610.3516;
sl=1/sf*1000;

 
load('mean_AER_all_trials.mat')
load('AER_reference_value_noPM.mat')
%%
time=[];
    for j= 1:size(mean_AER,2)
        time(j)=sl*j;
    end
    time=time-200;
 %%   
reftest=ones(1,122);    
semilogy(reftest, AER_max_peak2, 'ko') 
set(gca,'xtick',[])
box off
%%
how='AERpeaksusedforStandardization';
export_fig(how, '-bmp');
%% plot for AER on 1st tone


%dataAER=Alldat(:,1:489); %baseline + 600ms
%dataAER=Alldat(:,1099:1587); %last tone -200ms + 600ms

dataAER=((mean_AER.*100)./AER_max_peak2');
%%

 %time=time(1:489); %baseline + 600ms
%time=time(1099:1587); %last tone -200ms + 600ms

fill1=mean(dataAER,1)-std(dataAER);%/sqrt(length(data)); %for calculation of standard error
fill2=mean(dataAER,1)+std(dataAER);%/sqrt(length(data));

set(0,'defaultfigurecolor',[1 1 1])
figure
plot(time, mean(dataAER), 'color', grey ,'LineWidth', 3)
hold on
%plot(time, dataAER, 'k', 'LineWidth', 1)
% fill([time, flip(time)],[fill1';flip(fill2')], 'k', 'FaceAlpha','0.3','linestyle','none');
% hold on
%plot([0,0],[0,18e-15], 'k', 'LineWidth',2)
plot([0,0],[20 70], 'k', 'LineWidth',2)
hold on
plot([600,600],[20 70], 'k', 'LineWidth',2)
hold on
plot([1200,1200],[20 70], 'k', 'LineWidth',2)
hold on
plot([1800,1800],[20 70], 'k', 'LineWidth',2)
% plot([0,0],[-20 100], 'b', 'LineWidth',2)
% hold on
% plot([600,600],[-20 100], 'b', 'LineWidth',2)
% hold on
% plot([1200,1200],[-20 100], 'b', 'LineWidth',2)
% hold on
% plot([1800,1800],[-20 100], 'b', 'LineWidth',2)
xlabel('time (ms)')
ylabel('percentage change from first AER');
set(gca,'fontsize',14)
box off

figure
plot(time, dataAER(1,:))
hold on
plot([time(1,1), time(1,end)],[100,100])


%%
%all_selected_datasets=[selected_datasets1A; selected_datasets2A; selected_datasets1B; selected_datasets2B];


Info_1A=readtable('DatasetsIncludedInBrainAnalysis1Anew.csv');
Info_1B=readtable('DatasetsIncludedInBrainAnalysis1Bnew.csv');
Info_2A=readtable('DatasetsIncludedInBrainAnalysis2Anew.csv');
Info_2B=readtable('DatasetsIncludedInBrainAnalysis2Bnew.csv');
%% HRV parameters
HRVa=[Info_1A(:,8:17);Info_2A(:,8:17)];
HRVb=[Info_1B(:,8:17);Info_2B(:,8:17);];
HRVnames=HRVa.Properties.VariableNames;
%stim cond and sex
StimCona=[Info_1A(:,4);Info_2A(:,4)];
StimConb=[Info_1B(:,4);Info_2B(:,4);];

Sexa=[Info_1A(:,6);Info_2A(:,6)];
Sexb=[Info_1B(:,6);Info_2B(:,6);];

MeasurementNra=[Info_1A(:,7);Info_2A(:,7)];
MeasurementNrb=[Info_1B(:,7);Info_2B(:,7);];
% gestational age
GA_1A=table2array(Info_1A(:,5)); 
GA_1B=table2array(Info_1B(:,5)); 
GA_2A=table2array(Info_2A(:,5)); 
GA_2B=table2array(Info_2B(:,5)); 

GA_all=[GA_1A; GA_2A; GA_1B; GA_2B];% GA_1A; GA_2A; GA_1B; GA_2B];
% ID
ID_1A=table2array(Info_1A(:,2)); 
ID_1B=table2array(Info_1B(:,2)); 
ID_2A=table2array(Info_2A(:,2)); 
ID_2B=table2array(Info_2B(:,2)); 

ID_all=[ID_1A; ID_2A; ID_1B; ID_2B];
%%
ID_GA_all=[ID_all,GA_all];
[ID_sort, ID_idx]=sort(ID_all);
ID_GA_all_sort=ID_GA_all(ID_idx,:);
[unique_IDs, ia, ic]=unique(ID_GA_all_sort, 'rows');
%%

%%
load('wholetrial_standard.mat');
load('wholetrial_locglob.mat');
load('wholetrial_local.mat');
load('wholetrial_global.mat');

%% transform values to percentage values

AER_max_peakA=(AER_max_peak2(1,1:size(WTstandard,1)));
AER_max_peakB=(AER_max_peak2(1,size(WTstandard,1)+1:end));

dataGlob=((WTglobal.*100)./AER_max_peakB');
dataLoc=((WTlocal.*100)./AER_max_peakB');
dataLocglob=((WTlocglob.*100)./AER_max_peakA');
dataStd=((WTstandard.*100)./AER_max_peakA');
%% plot of AERs for tones 3 and 4
time_tone3=find(time>=1150&time<=1850);
time_tone4=find(time>=1750&time<=2450);
%time_tone4=find(time>=1850&time<=3000);
time_wind=find(time>=-50&time<=650);

tone3plotLg=dataLocglob(:,time_tone3(1,1:end-1));
tone4plotLg=dataLocglob(:,time_tone4);
diff_plotLg=mean(tone4plotLg-tone3plotLg);

tone3plotL=dataLoc(:,time_tone3(1,1:end-1));
tone4plotL=dataLoc(:,time_tone4);
diff_plotL=mean(tone4plotL-tone3plotL);

tone3plotG=dataGlob(:,time_tone3(1,1:end-1));
tone4plotG=dataGlob(:,time_tone4);
diff_plotG=mean(tone4plotG-tone3plotG);

tone3plotS=dataStd(:,time_tone3(1,1:end-1));
tone4plotS=dataStd(:,time_tone4);
diff_plotS=mean(tone4plotS-tone3plotS);

diff_plotGlob=(diff_plotG+diff_plotLg)/2;
diff_plotGlobStd=(diff_plotS+diff_plotL)/2;

diff_plotLoc=(diff_plotL+diff_plotLg)/2;
diff_plotLocStd=(diff_plotS+diff_plotG)/2;

time_plot=time(1,time_wind);
figure
plot(time_plot, diff_plotGlob, 'color', gr,  'LineWidth', 3)
hold on
plot(time_plot, diff_plotGlobStd, 'color', gs,  'LineWidth', 3)
hold on
% plot(time_plot, diff_plotG, 'color', g,  'LineWidth', 3)
% hold on
% plot(time_plot, diff_plotS, 'color', s,  'LineWidth', 3)
% hold on
plot([0,0],[-30 30], 'k', 'LineWidth',2)
xlabel('time (ms)')
xlim([-50 650]);
ylabel('difference T4-T3 in percent change');
set(gca,'fontsize',14)
legend('Global deviant', 'Global standard', 'Location', 'northoutside')
box off
%%
how='StandardT3T4';
export_fig(how, '-bmp');

%%
Globavg1=nanmean(dataGlob(:,AERtime1),2);
Globavg1sec=nanmean(dataGlob(:,AERtime1sec),2);
Globavg2=nanmean(dataGlob(:,AERtime2),2);
Globavg2sec=nanmean(dataGlob(:,AERtime2sec),2);
Globavg3=nanmean(dataGlob(:,AERtime3),2);
Globavg3sec=nanmean(dataGlob(:,AERtime3sec),2);
Globavg4=nanmean(dataGlob(:,AERtime4),2);
lateGlobavg1=nanmean(dataGlob(:,lateAERtime1),2);
lateGlobavg2=nanmean(dataGlob(:,lateAERtime2),2);
%lateGlobavg3=nanmean(dataGlob(:,lateAERtime3),2);
%%
figure
plot([mean(Globavg1), mean(Globavg1sec), mean(Globavg2), mean(Globavg2sec), mean(Globavg3), mean(Globavg3sec), mean(Globavg4), mean(lateGlobavg1), mean(lateGlobavg2)])
figure
boxplot([Globavg1, Globavg2, (Globavg3), (Globavg4), lateGlobavg1, lateGlobavg2, lateGlobavg3])
%%
Locavg1=nanmean(dataLoc(:,AERtime1),2);
Locavg1sec=nanmean(dataLoc(:,AERtime1sec),2);
Locavg2=nanmean(dataLoc(:,AERtime2),2);
Locavg2sec=nanmean(dataLoc(:,AERtime2sec),2);
Locavg3=nanmean(dataLoc(:,AERtime3),2);
Locavg3sec=nanmean(dataLoc(:,AERtime3sec),2);
Locavg4=nanmean(dataLoc(:,AERtime4),2);
lateLocavg1=nanmean(dataLoc(:,lateAERtime1),2);
lateLocavg2=nanmean(dataLoc(:,lateAERtime2),2);
%lateLocavg3=nanmean(dataLoc(:,lateAERtime3),2);
%%
figure
plot([mean(Locavg1), mean(Locavg1sec), mean(Locavg2), mean(Locavg2sec), mean(Locavg3), mean(Locavg3sec), mean(Locavg4), mean(lateLocavg1), mean(lateLocavg2)])
figure
boxplot([Locavg1, Locavg2, (Locavg3), (Locavg4), lateLocavg1, lateLocavg2, lateLocavg3])
%%
Locglobavg1=nanmean(dataLocglob(:,AERtime1),2);
Locglobavg1sec=nanmean(dataLocglob(:,AERtime1sec),2);
Locglobavg2=nanmean(dataLocglob(:,AERtime2),2);
Locglobavg2sec=nanmean(dataLocglob(:,AERtime2sec),2);
Locglobavg3=nanmean(dataLocglob(:,AERtime3),2);
Locglobavg3sec=nanmean(dataLocglob(:,AERtime3sec),2);
Locglobavg4=nanmean(dataLocglob(:,AERtime4),2);
lateLocglobavg1=nanmean(dataLocglob(:,lateAERtime1),2);
lateLocglobavg2=nanmean(dataLocglob(:,lateAERtime2),2);
%lateLocglobavg3=nanmean(dataLocglob(:,lateAERtime3),2);
%%
figure
plot([mean(Locglobavg1), mean(Locglobavg1sec), mean(Locglobavg2), mean(Locglobavg2sec), mean(Locglobavg3), mean(Locglobavg3sec), mean(Locglobavg4), mean(lateLocglobavg1), mean(lateLocglobavg2)])
figure
boxplot([Locglobavg1, Locglobavg2, (Locglobavg3), (Locglobavg4), lateLocglobavg1, lateLocglobavg2, lateLocglobavg3])
%%
Stdavg1=nanmean(dataStd(:,AERtime1),2);
Stdavg1sec=nanmean(dataStd(:,AERtime1sec),2);
Stdavg2=nanmean(dataStd(:,AERtime2),2);
Stdavg2sec=nanmean(dataStd(:,AERtime2sec),2);
Stdavg3=nanmean(dataStd(:,AERtime3),2);
Stdavg3sec=nanmean(dataStd(:,AERtime3sec),2);
Stdavg4=nanmean(dataStd(:,AERtime4),2);
lateStdavg1=nanmean(dataStd(:,lateAERtime1),2);
lateStdavg2=nanmean(dataStd(:,lateAERtime2),2);
%lateStdavg3=nanmean(dataStd(:,lateAERtime3),2);
%%
figure
plot([mean(Stdavg1), mean(Stdavg1sec), mean(Stdavg2), mean(Stdavg2sec), mean(Stdavg3), mean(Stdavg3sec), mean(Stdavg4), mean(lateStdavg1), mean(lateStdavg2)])
figure
boxplot([Stdavg1, Stdavg2, (Stdavg3), (Stdavg4), lateStdavg1, lateStdavg2, lateStdavg3])
%%
%% create table with subject characteristics
SubjChar=[repelem(unique_IDs(:,1),2),repelem(unique_IDs(:,2),2)] ;
%douplicate each unique ID and Ga, as theoretically every person has two
%recording blocks. 

BlockID=repmat(['A';'B'],size(unique_IDs,1),1); %create a vector with blocks A and B
%combine nfor in table and add a column of NaN's for eachv ariable that
%wants to be added
nancol=array2table(nan(size(SubjChar,1),49)); %depends on lize of added variables
SubjChar=[array2table(SubjChar), table(BlockID),nancol];
%name variables accordingly
SubjChar.Properties.VariableNames=[{'ID', 'GA','Block', 'T1_std', 'T1_std_sec','T2_std', 'T2_std_sec',...
    'T3_std', 'T3_std_sec','T4_std', 'T4_MM_std', 'T4_late_std', 'T1_loc', 'T1_loc_sec', 'T2_loc', 'T2_loc_sec',...
    'T3_loc', 'T3_loc_sec','T4_loc', 'T4_MM_loc', 'T4_late_loc', 'T1_glob', 'T1_glob_sec', 'T2_glob', 'T2_glob_sec',...
    'T3_glob', 'T3_glob_sec','T4_glob', 'T4_MM_glob', 'T4_late_glob', 'T1_locglob', 'T1_locglob_sec', 'T2_locglob', 'T2_locglob_sec',...
    'T3_locglob', 'T3_locglob_sec', 'T4_locglob', 'T4_MM_locglob', 'T4_late_locglob', 'StimCon', 'Sex', 'MeasurementNr'},  HRVnames];
%% A: combine values obtained in the analysis in the same fashion to match the subject characteristics table 
%create seperate table for each block type as number of recordings differ. 
blockA=table(repmat(['A'], [size(Info_1A,1)+size(Info_2A,1)], 1));
TabA=[[Info_1A(:,2); Info_2A(:,2)], ...
    [Info_1A(:,5); Info_2A(:,5)], blockA, array2table(Stdavg1), array2table(Stdavg1sec) ,array2table(Stdavg2),array2table(Stdavg2sec),...
    array2table(Stdavg3), array2table(Stdavg3sec), array2table(Stdavg4), array2table(lateStdavg1), array2table(lateStdavg2), ...
    array2table(Locglobavg1),array2table(Locglobavg1sec), array2table(Locglobavg2)...
    array2table(Locglobavg2sec),array2table(Locglobavg3), array2table(Locglobavg3sec), array2table(Locglobavg4), array2table(lateLocglobavg1), array2table(lateLocglobavg2)];
%name variables according to Subj charcteristics table
TabA.Properties.VariableNames={'ID', 'GA', 'Block', 'T1_std', 'T1_std_sec','T2_std', 'T2_std_sec',...
    'T3_std', 'T3_std_sec','T4_std', 'T4_MM_std', 'T4_late_std', 'T1_locglob', 'T1_locglob_sec', 'T2_locglob', 'T2_locglob_sec',...
    'T3_locglob', 'T3_locglob_sec', 'T4_locglob', 'T4_MM_locglob', 'T4_late_locglob'};
TabA=[TabA, StimCona, Sexa, MeasurementNra, HRVa];
%% B: combine values obtained in the analysis in the same fashion to match the subject characteristics table 
%create seperate table for each block type as number of recordings differ. 
blockB=table(repmat(['B'], [size(Info_1B,1)+size(Info_2B,1)], 1));
TabB=[[ Info_1B(:,2); Info_2B(:,2)], ...
    [Info_1B(:,5); Info_2B(:,5)], blockB, array2table(Locavg1), array2table(Locavg1sec), array2table(Locavg2), array2table(Locavg2sec),...
    array2table(Locavg3), array2table(Locavg3sec), array2table(Locavg4), array2table(lateLocavg1), ...
    array2table(lateLocavg2), array2table(Globavg1), array2table(Globavg1sec), array2table(Globavg2), ...
    array2table(Globavg2sec),array2table(Globavg3), array2table(Globavg3sec), array2table(Globavg4), array2table(lateGlobavg1), array2table(lateGlobavg2)];
%name variables according to Subj charcteristics table
TabB.Properties.VariableNames={'ID', 'GA', 'Block','T1_loc', 'T1_loc_sec', 'T2_loc', 'T2_loc_sec',...
    'T3_loc', 'T3_loc_sec','T4_loc', 'T4_MM_loc', 'T4_late_loc', 'T1_glob', 'T1_glob_sec', 'T2_glob', 'T2_glob_sec',...
    'T3_glob', 'T3_glob_sec','T4_glob', 'T4_MM_glob', 'T4_late_glob'};
TabB=[TabB, StimConb, Sexb, MeasurementNrb, HRVb];
%% fill in values from Table A into SubjChar
for kj=1:size(TabA,1) %loop over the size of the table with the MEG data
    for tj=1:size(SubjChar,1) %second lop to check each entry of the Subj char table
kidx1(tj,:)=SubjChar{tj,1}== TabA{kj, 1}; %check for matching IDs
kidx2(tj,:)=SubjChar{tj,2}== TabA{kj, 2}; %matching GA's
kidx3(tj,:)=SubjChar{tj,3}== TabA{kj, 3}; %matching block types
    end
    kidx=sum([kidx1, kidx2, kidx3],2); %sum matches up
    ridx=find(kidx==3); %thre should be only one occation, where all three match
    SubjChar.T1_std(ridx)=TabA.T1_std(kj); %take that spot and fill in the value from the MEG table
     SubjChar.T1_std_sec(ridx)=TabA.T1_std_sec(kj);
    SubjChar.T2_std(ridx)=TabA.T2_std(kj);
    SubjChar.T2_std_sec(ridx)=TabA.T2_std_sec(kj);
    SubjChar.T3_std(ridx)=TabA.T3_std(kj);
    SubjChar.T3_std_sec(ridx)=TabA.T3_std_sec(kj);
    SubjChar.T4_std(ridx)=TabA.T4_std(kj);
    SubjChar.T4_MM_std(ridx)=TabA.T4_MM_std(kj);
    SubjChar.T4_late_std(ridx)=TabA.T4_late_std(kj);
    SubjChar.T1_locglob(ridx)=TabA.T1_locglob(kj); %take that spot and fill in the value from the MEG table
    SubjChar.T1_locglob_sec(ridx)=TabA.T1_locglob_sec(kj);
    SubjChar.T2_locglob(ridx)=TabA.T2_locglob(kj);
    SubjChar.T2_locglob_sec(ridx)=TabA.T2_locglob_sec(kj);
    SubjChar.T3_locglob(ridx)=TabA.T3_locglob(kj);
    SubjChar.T3_locglob_sec(ridx)=TabA.T3_locglob_sec(kj);
    SubjChar.T4_locglob(ridx)=TabA.T4_locglob(kj);
    SubjChar.T4_MM_locglob(ridx)=TabA.T4_MM_locglob(kj);
    SubjChar.T4_late_locglob(ridx)=TabA.T4_late_locglob(kj);
    SubjChar.StimCon(ridx)=TabA.StimCondition(kj);
    SubjChar.Sex(ridx)=char(TabA.gender(kj));
    SubjChar(ridx, 42:52)=TabA(kj,24:34);
end
%%
%% fill in values from Table B into SubjChar
for kj=1:size(TabB,1) %loop over the size of the table with the MEG data
    for tj=1:size(SubjChar,1) %second lop to check each entry of the Subj char table
kidx1(tj,:)=SubjChar{tj,1}== TabB{kj, 1}; %check for matching IDs
kidx2(tj,:)=SubjChar{tj,2}== TabB{kj, 2}; %matching GA's
kidx3(tj,:)=SubjChar{tj,3}== TabB{kj, 3}; %matching block types
    end
    kidx=sum([kidx1, kidx2, kidx3],2); %sum matches up
    ridx=find(kidx==3); %thre should be only one occation, where all three match
    SubjChar.T1_loc(ridx)=TabB.T1_loc(kj); %take that spot and fill in the value from the MEG table
    SubjChar.T1_loc_sec(ridx)=TabB.T1_loc_sec(kj);
    SubjChar.T2_loc(ridx)=TabB.T2_loc(kj);
    SubjChar.T2_loc_sec(ridx)=TabB.T2_loc_sec(kj);
    SubjChar.T3_loc(ridx)=TabB.T3_loc(kj);
    SubjChar.T3_loc_sec(ridx)=TabB.T3_loc_sec(kj);
    SubjChar.T4_loc(ridx)=TabB.T4_loc(kj);
    SubjChar.T4_MM_loc(ridx)=TabB.T4_MM_loc(kj);
    SubjChar.T4_late_loc(ridx)=TabB.T4_late_loc(kj);
    SubjChar.T1_glob(ridx)=TabB.T1_glob(kj); %take that spot and fill in the value from the MEG table
    SubjChar.T1_glob_sec(ridx)=TabB.T1_glob_sec(kj);
    SubjChar.T2_glob(ridx)=TabB.T2_glob(kj);
    SubjChar.T2_glob_sec(ridx)=TabB.T2_glob_sec(kj);
    SubjChar.T3_glob(ridx)=TabB.T3_glob(kj);
    SubjChar.T3_glob_sec(ridx)=TabB.T3_glob_sec(kj);
    SubjChar.T4_glob(ridx)=TabB.T4_glob(kj);
    SubjChar.T4_MM_glob(ridx)=TabB.T4_MM_glob(kj);
    SubjChar.T4_late_glob(ridx)=TabB.T4_late_glob(kj);
     SubjChar.StimCon(ridx)=TabB.StimCondition(kj);
    SubjChar.Sex(ridx)=char(TabB.gender(kj));
    SubjChar(ridx, 42:52)=TabB(kj,24:34);
end
%%
writetable(SubjChar, 'AVG_indConditions_pct_inclsec_FAIRY.csv' )
%%
%%
 %Data Export for sharing
 % all trials
 Alltrials=num2cell(mean_AER,2);
 %maximum peaks
 Peakmax=AER_max_peak2';
 % all trials normalized
 Alltrials_normalized=num2cell(dataAER,2);
 %non normalized data
cellglobal=num2cell([WTlocglob;WTglobal], 2); 
cellloc=num2cell([WTstandard; WTlocal],2);
%normlaized data
cellglobal_norm=num2cell([dataLocglob;dataGlob], 2); 
cellloc_norm=num2cell([dataStd; dataLoc],2);

Blocklabel=[repmat('A (ssss)',size(WTstandard,1), 1); repmat('B (sssd)',size(WTlocal,1), 1)];

Datatable=table(ID_all, GA_all, Blocklabel, ...
    Alltrials, Peakmax, Alltrials_normalized, ...
    cellloc, cellloc_norm, cellglobal, cellglobal_norm);
Datatable.Properties.VariableNames={'ID', 'GA', 'Blocklabel', ...
    'all_trials', 'maximum_peak_amplitude', 'all_trails_normalized', 'global_standards', ...
    'global_standards_normalized', 'global_deviants', 'global_deviants_normalized'};

save('Data_traces.mat', 'Datatable')