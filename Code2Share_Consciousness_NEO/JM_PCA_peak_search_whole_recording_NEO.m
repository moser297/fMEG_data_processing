function [high5, high52, high53, high5h]=JM_PCA_peak_search_whole_recording_NEO(data_new_filtered, data_new_unfiltered, path, channelinfo, sf, filename)
% this function is used to create plos that help to evaluate data quality and especially
%amount and location of leftover heart activity. To judge leftover hesrt
%activity, markers from heart removal are used and averaged
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

%% variance over whole recording

colorA=var(data_new_filtered');


%color=i_coeff1(4,:);
figure
scatter(cell2mat(channelinfo(:,2)), cell2mat(channelinfo(:,3)),200, colorA(1,:), 'filled')
title('Variance on each channel over whole recording')
colorbar
saveas(figure(1), [filename_short, '_var'], 'bmp')

%% Average leftover heart signal
addpath(genpath('/opt_prg/fieldtrip-20170202'), '-end') %add path to the end of the search path as 
%fieldtrip also has a pca function, but we want to use the function from the stats toolbox

cd(path)
event = ft_read_event(filename);
%make sure, teht event vector has intended orientation as orientations
%differ between FAIRY processing and voncentional processing
if size(event,1)<size(event,2)
    event=event';
end
%average over fMCG trigger for getting left over heart signal
%NEO: only fMCG!
triggernameF='FLORA_fMCG';
triggername1='fMCG2';

triggername3='fMCG';


triggervec1 = [];
samples = [];
%check, if there is a FLORA trigger
for i = 2:size(event,1)
    type = event(i).type;
    samples = [samples event(i).sample];
    if strcmp(type, triggernameF)
        triggervec1 = [triggervec1 event(i).sample];
    end
end
%if not, is there a fMCG2 trigger?
if isempty(triggervec1)
    for i = 2:size(event,1)
        type = event(i).type;
        samples = [samples event(i).sample];
            if strcmp(type, triggername1)
                triggervec1 = [triggervec1 event(i).sample];
            end
    end
end
%check, if there are actually fMCG2 trigger in the dataset. Otherwise use
%fMCG
if isempty(triggervec1)
    for i = 2:size(event,1)
        type = event(i).type;
        samples = [samples event(i).sample];
            if strcmp(type, triggername3)
                triggervec1 = [triggervec1 event(i).sample];
            end
    end
end
triggervec1 = unique(triggervec1);
%cut triggervec to 6 min length of dataset
[~, indy]=find(triggervec1 >= size(data_new_unfiltered,2));
if ~isempty(indy)
triggervec1=triggervec1(:,1:indy-1);
end
%triggervec1=fetal heart leftover


sl=1/sf*1000; %length of 1 sample 

t1=round(100/sl); %100ms before heartbeat in ms
t2=round(200/sl); %200ms after heart beat in ms

%leftover fetal heart
for h=2:size(triggervec1,2)-2 %not taking first and last beat into account as there is possibly not enough time around them
    leftover_fetal_h(:,:,h)=data_new_unfiltered(:,triggervec1(1,h)-t1:triggervec1(1,h)+t2);
end

%% baselinecorrection 
for k=1:size(leftover_fetal_h,3)
    base=mean(leftover_fetal_h(:,1:t1,k),2);
    leftover_fetal_h_b(:,:,k)=leftover_fetal_h(:,:,k)-base;
end

mean_leftover_fetal_h_b=mean(leftover_fetal_h_b, 3);

%% plot of leftover heart
figure
for ij=1:size(mean_leftover_fetal_h_b, 1)
    plot(mean_leftover_fetal_h_b(ij,:))
    hold on
    title('Fetal heart leftover')
end
saveas(figure(2), [filename_short, '_Fetal_heart_leftover'], 'bmp')

%% select channels with highest amplitude and plot them


%color1=max(mean_leftover_fetal_h_b'); %max amplitude sorted by channel
% amplitude around R peak - change 30.07.2019
color1=max(abs(mean_leftover_fetal_h_b(:,47:77)')); % values represent 25ms before and 25ms after R-peak timepoint (timepoint 62)
[~, indexh] = sort(color1, 'descend');
% 5 channels most important for PC1
high5h=indexh(1,1:5);

figure
set(gcf, 'Position', [1300, 450, 1300, 450])
subplot(1,2,1)
scatter(cell2mat(channelinfo(:,2)), cell2mat(channelinfo(:,3)),200, color1(1,:), 'filled')
title('channels with most fetal heart leftover')
colorbar
subplot(1,2,2)
scatter(cell2mat(channelinfo(:,2)), cell2mat(channelinfo(:,3)),200, color1(1,:), 'filled')
title('5 channels with most fetal heart leftover')
colorbar
hold on 
scatter(cell2mat(channelinfo(high5h,2)),cell2mat(channelinfo(high5h,3)),200, 'filled', 'r')
saveas(figure(3), [filename_short, '_channels_fetal_heart_leftover'], 'bmp')

%% PCA
%score=signal over time of each principal component
%coeff=how much influence does each cahnnel have on the each component
%latent how much variance does each component explain (how important is each component)
[coeff,score,~, ~, explained] = pca(data_new_filtered');
%inverse matrix 
i_coeff=inv(coeff);
%inverse matrix for first 4 components -i_coeff_sel: how much influence
%does the component have on each channel
i_coeff_sel=i_coeff(1:4,:);
%selected data is data that only consists of first components, everything
%else is filtered out
selected_data = i_coeff_sel(1,:)'* score(:,1)';
selected_data_2 = i_coeff_sel(2,:)'* score(:,2)';
selected_data_3 = i_coeff_sel(3,:)'* score(:,3)';
%% channel plot
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
set(0,'defaultfigurecolor',[1 1 1])
scatter(cell2mat(channelinfo(:,2)), cell2mat(channelinfo(:,3)),200, color(1,:), 'filled')
title(['whole recording PC1, Explained Variance:' num2str(explained(1,1))])
colorbar
subplot(1,2,2)
scatter(cell2mat(channelinfo(:,2)), cell2mat(channelinfo(:,3)),200, color(1,:), 'filled')
title('5 most important channels for PC1')
colorbar
hold on 
scatter(cell2mat(channelinfo(high5,2)),cell2mat(channelinfo(high5,3)),200, 'filled', 'r')
saveas(figure(4), [filename_short, '_PC1_whole'], 'bmp')

% PC2
color2=i_coeff(2,:);

coeff_PC2=coeff(:,2);
%sort coefficients of PC1
[~, index2] = sort(coeff_PC2, 'descend');
% 5 channels most important for PC2
high52=index2(1:5,1);

figure
set(gcf, 'Position', [1300, 450, 1300, 450])
subplot(1,2,1)
scatter(cell2mat(channelinfo(:,2)), cell2mat(channelinfo(:,3)),200, color2(1,:), 'filled')
title(['whole recording PC2, Explained Variance:' num2str(explained(2,1))])
colorbar
subplot(1,2,2)
scatter(cell2mat(channelinfo(:,2)), cell2mat(channelinfo(:,3)),200, color2(1,:), 'filled')
title('5 most important channels for PC2')
colorbar
hold on 
scatter(cell2mat(channelinfo(high52,2)),cell2mat(channelinfo(high52,3)),200, 'filled', 'r')
saveas(figure(5), [filename_short, '_PC2_whole'], 'bmp')

% PC3
color3=i_coeff(3,:);

coeff_PC3=coeff(:,3);
%sort coefficients of PC1
[~, index3] = sort(coeff_PC3, 'descend');
% 5 channels most important for PC2
high53=index3(1:5,1);

figure
set(gcf, 'Position', [1300, 450, 1300, 450])
subplot(1,2,1)
scatter(cell2mat(channelinfo(:,2)), cell2mat(channelinfo(:,3)),200, color3(1,:), 'filled')
title(['whole recording PC3, Explained Variance:' num2str(explained(3,1))])
colorbar
subplot(1,2,2)
scatter(cell2mat(channelinfo(:,2)), cell2mat(channelinfo(:,3)),200, color3(1,:), 'filled')
title('5 most important channels for PC3')
colorbar
hold on 
scatter(cell2mat(channelinfo(high53,2)),cell2mat(channelinfo(high53,3)),200, 'filled', 'r')
saveas(figure(6), [filename_short, '_PC3_whole'], 'bmp')
%% signal traces
sl=1000/sf;
time=sl;
for i=1:4999
time(i+1)=time(i) +sl;
end
figure
set(gcf, 'Position', [1000, 550, 1000, 550])
subplot(3,1,1)
plot(time,selected_data(:,1001:6000)') %random time window from data
title('PC1 trace whole')
subplot(3,1,2)
plot(time,selected_data_2(:,1001:6000)')
title('PC2 trace whole')
subplot(3,1,3)
plot(time,selected_data_3(:,1001:6000)')
title('PC3 trace whole')
saveas(figure(7), [filename_short, '_trace_PCA'], 'bmp')

%%
close all

end