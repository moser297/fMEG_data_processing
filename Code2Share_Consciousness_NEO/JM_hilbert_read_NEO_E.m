function [channelinfo,data_new_filtered, data_new_heart_dc, sf]=JM_hilbert_read_NEO_E(path, filename)
%This function reads in data and gives out filtered data for processing
%(1-10Hz), for averaging of heart leftover, data is not filtered but the
%offset is removed
%path is the path to the folder of the subject data
%filename is the name of the data file processed with Hilbert/DataEditor
% for comparability, data is reduced to first 6 minutes of recording

%%
cd(path)
%readCTFdas reads in all information about a data file
ds=readCTFds(filename);
% ds_heart=readCTFds(filename);
%getCTFdata: 1, data 2. [] read complete dataset 3. [] read all channels 4. use
%unit Tesla - this function is needed to access the actual MEG data. 
d1=getCTFdata(ds,[],[],'t');
%chec, if data consists of several trials
if numel(size(d1))==3
    d2=permute(d1, [1 3 2]); %data needs to be reshaped as it consists of two trials, first dimensions are permuted to achive the correct result for reshaping
    d3=reshape(d2, [size(d1,1)*2, size(d1,2)]);
else 
    d3=d1;
end
    
% d_heart1=getCTFdata(ds_heart,[],[],'t');
% d_heart2=permute(d_heart1, [1 3 2]);
% d_heart3=reshape(d_heart2, [size(d_heart1,1)*2, size(d_heart1,2)]);
%% shortening Data to 6 minutes for comparability with OGTT data - not applied anymore

% sp=round(ds.res4.sample_rate*360); %360 seconds for 6min time window
% if size(d3,1)>sp
%     d=d3(1:sp,:);
% else
    d=d3;
% end
%% channel position 
%extract info about sensor positions from ds structure
sensorinfo=ds.res4.senres;
%get numbers of all MEG cahnnels
MEGchan=startsWith(ds.res4.chanNames, 'M');
I=find(MEGchan, 1, 'first'); %first MEG channel
J=find(MEGchan, 1, 'last'); %last MEG channel

% get position of every MEG channel (out=outer coil of gradiometer, in=
% inner coil of gradiometer). build mean position
for i=I:J
    pos0_out=sensorinfo(i).pos0(:,1)';
    pos0_in=sensorinfo(i).pos0(:,2)';
    pos0_mean(i,:)=mean([pos0_out; pos0_in],1);
end
%delete 0 rows from mean
pos0_mean(all(~pos0_mean,2),:)=[];
%only MEG channels
d=d(:,I:J);
%% in this version all channels are used and not only the relevant ones from the previous analysis
%d_heart=d_heart(:,I:J);
% relevant channel positions same as FAIRY: 
%relevantChannels=[3;4;5;6;7;8;9;10;11;21;24;25;26;27;28;29;30;31;32;33;34;35;36;37;38;39;40;41;42;43;44;45;46;47;48;49;50;51;52;53;54;55;56;58;59;60;61;91;94;95;96;97;98;99;100;101;102;103;104;105;106;107;108;109;110;111;112;113;114;115;116;117;118;119;120;121;122;123;124;125;126;128;129;130;131];
%reduction of channels to 64 channels as all subjects are in SL
%relevantChannels=[3;4;5;6;7;8;9;10;21;24;25;26;27;28;29;30;31;32;33;34;35;36;37;38;39;40;41;42;43;44;45;46;47;48;49;50;91;94;95;96;97;98;99;100;101;102;103;104;105;106;107;108;109;110;111;112;113;114;115;116;117;118;119;120];
% names of 64 relevant channels to correctly account for missing channels
% relevantChannelNames=['MCC0'; 'MCD0'; 'MCE0'; 'MCF0'; 'MCG0'; 'MCH0'; 'MCI0'; 'MCJ0'; 'MLC1'; 'MLD1'; 'MLD2'; 'MLD3'; 'MLE1';...
% 'MLE2'; 'MLE3'; 'MLF1'; 'MLF2'; 'MLF3'; 'MLG1'; 'MLG2'; 'MLG3'; 'MLH1'; 'MLH2'; 'MLH3'; 'MLH4'; 'MLI1';...
% 'MLI2'; 'MLI3'; 'MLI4'; 'MLI5'; 'MLJ1'; 'MLJ2'; 'MLJ3'; 'MLJ4'; 'MLJ5'; 'MLJ6'; 'MRC1';...
% 'MRD1'; 'MRD2'; 'MRD3'; 'MRE1'; 'MRE2'; 'MRE3'; 'MRF1'; 'MRF2'; 'MRF3'; 'MRG1'; 'MRG2'; 'MRG3'; 'MRH1';...
% 'MRH2'; 'MRH3'; 'MRH4'; 'MRI1'; 'MRI2'; 'MRI3'; 'MRI4'; 'MRI5'; 'MRJ1'; 'MRJ2'; 'MRJ3'; 'MRJ4'; 'MRJ5'; 'MRJ6'];                       
% 
names_all_chan=ds.res4.chanNames;
names_MEG_chan=names_all_chan(MEGchan,:);
% 
% empty_list=zeros(size(names_MEG_chan,1),1);
% ind=[];
% for i=1:size(relevantChannelNames, 1)
%     present_chan=contains(names_MEG_chan, relevantChannelNames(i,:));
%     if i==1
%      relevantChannels=empty_list+present_chan;
%     else 
%         relevantChannels=relevantChannels+present_chan;
%     end
%     if present_chan==0
%         ind(i)=i;
%     end
% end
% 
% relevantChannels=logical(relevantChannels);
% 
% channelinfo_pos=pos0_mean(relevantChannels,:);
% data_new=d(:,relevantChannels);
% data_new=data_new';
channelinfo_pos=pos0_mean;
data_new=d';
data_new_heart=data_new; %extra dataset for heart as it is not bandpassfiltered and offset is removed
% data_new_heart=d_heart(:,relevantChannels);
% data_new_heart=data_new_heart';
%% removal of offset for heart data
base=mean(data_new_heart,2);
data_new_heart_dc=data_new_heart-base;

%% filter 
%format of data_new: number of channels x data length!
%NEO data: 1-15Hz filter
sf=ds.res4.sample_rate;
data_new_filtered=KS_butterworthfilter(data_new, sf, 1 ,15);

%save channel names with channel info
%1. detect relevant channels that were not present in dataset
% ind(ind==0)=[];
% presentChannelNames=cellstr(relevantChannelNames);
% presentChannelNames(ind,:)=[];

%version including all channels:
presentChannelNames=cellstr(names_MEG_chan(:,1:4));

%2. combine names with cahnnel position
channelinfo=[presentChannelNames,num2cell(channelinfo_pos)];
%if channel check wants to be done
%is_ok = KS_check_channels(data_new_filtered,channelinfo_ext,5,path,1);
%%
end