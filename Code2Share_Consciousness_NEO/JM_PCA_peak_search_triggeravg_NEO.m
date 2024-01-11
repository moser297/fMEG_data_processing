function [whole_trial_m, remaining_trials, high5_m, high5_m2, high5_m3, explained, ...
    high5, high52, high53, pks, ptm, pks2, ptm2, pks3, ptm3, ...
     pks_nn, ptm_nn, pks2_nn, ptm2_nn, pks3_nn, ptm3_nn]=JM_PCA_peak_search_triggeravg_NEO(data_new_filtered, channelinfo, sf, path, filename, trig, tbeft, taft)
%this function is used to evaluate the sequences SbefD or Dev. -only 1
%trigger is averaged. In comparison to v2, the plusminus noise criterion is
%used.
 %path is the path to the folder of the subject data
%filename is the name of the original data file
%trig can be either DevA, DevB, DevC, DevD or Sbef and codes for the
% triger that wants to be evaluated
%time before trigger in ms: tbeft
%time after trigger in ms: taft
if isempty(strfind(filename,'fMEG_data.ds'))
filename_short=filename(1,1:18);
else
    filename_short=path(1,68:87);
end

set(0,'defaultfigurecolor',[1 1 1])
%% search for trigger events
addpath(genpath('/opt_prg/fieldtrip-20170202'), '-end') %add path to the end of the search path as 
%fieldtrip also has a pca function, but we want to use the function from the stats toolbox

cd(path)
event = ft_read_event(filename);
%make sure, teht event vector has intended orientation as orientations
%differ between FAIRY processing and voncentional processing
if size(event,1)<size(event,2)
    event=event';
end

triggervec = [];
samples = [];

if trig =='DevA'
triggername= 'a_500_500_500_500';
elseif trig =='DevB'
triggername= 'b_500_500_500_750';
elseif trig == 'DevC'
triggername= 'c_750_750_750_750';
elseif trig == 'DevD'
triggername= 'd_750_750_750_500';
elseif trig == 'Sbef'
triggername= 'SbefD';
end
%types = event(:).type;
for i = 2:size(event)
    type = event(i).type;
    samples = [samples event(i).sample];
    if strcmp(type, triggername)
        triggervec = [triggervec event(i).sample];
    end
end
triggervec = unique(triggervec);
% all 46 trials cut at the trigger. t1 samples before trigger onset and t2,
% samples after onset of tone; from triger to actual sound output there is a delay of 27ms

sl=1/sf*1000;
t1=round(tbeft/sl);
t2=round(taft/sl);
%in case data is cut shorter previously
[~, indy]=find(triggervec >= size(data_new_filtered,2)-t2);
if ~isempty(indy)
triggervec=triggervec(:,1:indy-1);
end

t_delay=round(27/sl); %delay of sound output from trigger to ear

triggervec_d=triggervec + t_delay; %trigger including delay of sound

%in case, first trigger is too early, so no baseline would be possible,
%remove trigger
if triggervec_d(1,1)<t1
    triggervec_d(:,1)=[];
end 

for h=1:size(triggervec,2)
    whole_trial(:,:,h)=data_new_filtered(:,triggervec_d(1,h)-t1:triggervec_d(1,h)+t2);
end
%% thresholddetect
kk=1;
while kk<= size(whole_trial,3)
    current=whole_trial(:,:,kk);
    if any(any(abs(current) >=1e-12))==1
    whole_trial(:,:,kk)=[];
    else 
        kk=kk+1;
    end
end
remaining_trials=size(whole_trial,3);
%% baselinecorrection 
for k=1:size(whole_trial,3)
    base=mean(whole_trial(:,1:t1,k),2);
    whole_trial_b(:,:,k)=whole_trial(:,:,k)-base;
end

%%
%mean over all trials
whole_trial_m=mean(whole_trial_b,3);
%% PCA
%score=signal over time of each principal component
%coeff=how much influence does each cahnnel have on the each component
%latent how much variance does each component explain (how important is each component)
[coeff,score,latent, tsquared, explained] = pca(whole_trial_m');
%inverse matrix 
i_coeff=inv(coeff);
%inverse matrix for first 4 components -i_coeff_sel: how much influence
%does the component have on each channel
i_coeff_sel=i_coeff(1:4,:);
%selected data is data that only consists of first components, everything
%else is filtered out
%PC1:
selected_data = i_coeff_sel(1,:)'* score(:,1)';
%PC2:
selected_data2 = i_coeff_sel(2,:)'* score(:,2)';
%PC3:
selected_data3 = i_coeff_sel(3,:)'* score(:,3)';
%% Find most important channels for PC1 and create channel plot
%plot most important channels for PC 1
color=i_coeff(1,:);

coeff_PC1=coeff(:,1);
%sort coefficients of PC1
[~, index] = sort(coeff_PC1, 'descend');
% 5 channels most important for PC1
high5=index(1:5,1);

figure
set(gcf, 'Position', [1300, 450, 1300, 450])
subplot(1,2,1)
scatter(cell2mat(channelinfo(:,2)), cell2mat(channelinfo(:,3)),200, color(1,:), 'filled')
title(['Trigger avg channels for PC1, Explained Variance:' num2str(explained(1,1))])
colorbar
subplot(1,2,2)
scatter(cell2mat(channelinfo(:,2)), cell2mat(channelinfo(:,3)),200, color(1,:), 'filled')
title('5 most important channels for PC1')
colorbar
hold on 
scatter(cell2mat(channelinfo(high5,2)),cell2mat(channelinfo(high5,3)),200, 'filled', 'r')
saveas(figure(1), [filename_short, trig, '_channelplot_PC1'], 'bmp')
 %% same for PC2

 color2=i_coeff(2,:);
 
coeff_PC2=coeff(:,2);
%sort coefficients of PC1
[~, index2] = sort(coeff_PC2, 'descend');
% 5 channels most important for PC2
high52=index2(1:5,1);
% 
figure
set(gcf, 'Position', [1300, 450, 1300, 450])
subplot(1,2,1)
scatter(cell2mat(channelinfo(:,2)), cell2mat(channelinfo(:,3)),200, color2(1,:), 'filled')
title(['Trigger avg channels for PC2, Explained Variance:' num2str(explained(2,1))])
colorbar
subplot(1,2,2)
scatter(cell2mat(channelinfo(:,2)), cell2mat(channelinfo(:,3)),200, color2(1,:), 'filled')
title('5 most important channels for PC2')
colorbar
hold on 
scatter(cell2mat(channelinfo(high52,2)),cell2mat(channelinfo(high52,3)),200, 'filled', 'r')
saveas(figure(2), [filename_short, trig, '_channelplot_PC2'], 'bmp')

 %% same for PC3

 color3=i_coeff(3,:);
 
coeff_PC3=coeff(:,3);
%sort coefficients of PC1
[~, index3] = sort(coeff_PC3, 'descend');
% 5 channels most important for PC2
high53=index3(1:5,1);
% 
figure
set(gcf, 'Position', [1300, 450, 1300, 450])
subplot(1,2,1)
scatter(cell2mat(channelinfo(:,2)), cell2mat(channelinfo(:,3)),200, color3(1,:), 'filled')
title(['Trigger avg channels for PC3, Explained Variance:' num2str(explained(3,1))])
colorbar
subplot(1,2,2)
scatter(cell2mat(channelinfo(:,2)), cell2mat(channelinfo(:,3)),200, color3(1,:), 'filled')
title('5 most important channels for PC3')
colorbar
hold on 
scatter(cell2mat(channelinfo(high53,2)),cell2mat(channelinfo(high53,3)),200, 'filled', 'r')
saveas(figure(3), [filename_short, trig, '_channelplot_PC3'], 'bmp')
%% plot channels over trial
% time over whole trial
for j= 1:size(whole_trial_m,2)
    time(j)=((1/sf)*1000)*j-tbeft;
end
% 
figure
set(gcf, 'Position', [700, 850, 700, 850])
subplot(6,1,1)
plot(time, selected_data(high5,:)')
title('Data weighted with PC1 - highest channels')
subplot(6,1,2)
plot(time, whole_trial_m(high5,:)')
title('Same channels (PC1) original data')
subplot(6,1,3)
plot(time, selected_data2(high52,:)')
title('Data weighted with PC2 - highest channels')
subplot(6,1,4)
plot(time, whole_trial_m(high52,:)')
title('Same channels (PC2) original data')
subplot(6,1,5)
plot(time, selected_data3(high53,:)')
title('Data weighted with PC3 - highest channels')
subplot(6,1,6)
plot(time, whole_trial_m(high53,:)')
title('Same channels (PC3) original data')
saveas(figure(4), [filename_short, trig, '_traces_PC1_PC2_PC3'], 'bmp')

%% peak search
%peak search on those channels, original data
%PC1:
%peak search on original data
high5_m=rms(whole_trial_m(high5,:));
%peak search on weighted data
%high5_m=rms(selected_data(high5,:));
[pks, ptm]=findpeaks(high5_m, 'MinPeakHeight', 5e-15);
%PC2:
%peak search on original data
high5_m2=rms(whole_trial_m(high52,:));
%peak search on weighted data
%high5_m=rms(selected_data(high5,:));
[pks2, ptm2]=findpeaks(high5_m2, 'MinPeakHeight', 5e-15);

%PC3:
%peak search on original data
high5_m3=rms(whole_trial_m(high53,:));
%peak search on weighted data
%high5_m=rms(selected_data(high5,:));
[pks3, ptm3]=findpeaks(high5_m3, 'MinPeakHeight', 5e-15);

%% peak search with +/- noise criterion
%% plusminus as noise criterion
%mean over all trials, if every second trial is subtracted
%number of itterations (a) is defined. In case there is an uneven number of
%trials, the last trial is not used. 
if mod(size(whole_trial_b,3),2)==0
    a=size(whole_trial_b,3);
else
    a=size(whole_trial_b,3)-1;
end

for k=1:a
    if mod(k,2)==0
    whole_trial_pm(:,:,k)=whole_trial_b(:,:,k).*-1;
    else
    whole_trial_pm(:,:,k)=whole_trial_b(:,:,k);    
    end
end
%mean over all trials
whole_trial_pm_m=mean(whole_trial_pm,3);
% rms over 5 important channels
%PC1:
high5_m_pm=rms(whole_trial_pm_m(high5,:));
%PC2:
high5_m2_pm=rms(whole_trial_pm_m(high52,:));
%PC3:
high5_m3_pm=rms(whole_trial_pm_m(high53,:));
%% comparison of plusminus with peaks PC1
pks_pm=high5_m_pm(ptm);
pks_nn=pks; %_nn -> no noise peaks left over
ptm_nn=ptm; 
ct=1;
while ct<=size(ptm_nn,2)
 if pks_pm(1,ct)>pks_nn(1,ct)/2 %if plus/minus peak is higher then half of the actual peak size - remove this peak
     pks_nn(1,ct)=NaN;
     ptm_nn(1,ct)=NaN;
 else
     ct=ct+1;
 end
end
pks_nn(isnan(pks_nn))=[];
ptm_nn(isnan(ptm_nn))=[];
%% comparison of plusminus with peaks PC2
pks_pm2=high5_m2_pm(ptm2);
pks2_nn=pks2;
ptm2_nn=ptm2; 
ct=1;
while ct<=size(ptm2_nn,2)
 if pks_pm2(1,ct)>pks2_nn(1,ct)/2
     pks2_nn(1,ct)=NaN;
     ptm2_nn(1,ct)=NaN;
 else
     ct=ct+1;
 end
end
pks2_nn(isnan(pks2_nn))=[];
ptm2_nn(isnan(ptm2_nn))=[];
%% comparison of plusminus with peaks PC3
pks_pm3=high5_m3_pm(ptm3);
pks3_nn=pks3;
ptm3_nn=ptm3; 
ct=1;
while ct<=size(ptm3_nn,2)
 if pks_pm3(1,ct)>pks3_nn(1,ct)/2
     pks3_nn(1,ct)=NaN;
     ptm3_nn(1,ct)=NaN;
 else
     ct=ct+1;
 end
end
pks3_nn(isnan(pks3_nn))=[];
ptm3_nn(isnan(ptm3_nn))=[];
%%
close all
end