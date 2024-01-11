function [MMbetween]=MMbetweenSequence_v3(path, sf, tbeft, resultsSbefD, resultsDev, PCA_check_SbefD, PCA_check_Dev)
%this function is similar to MMbetweenSequence_v2 but it takes the rms of
%the 5 selected channels instead of the mean


rms_PC1_SbefD=cell2mat(resultsSbefD{6,2}); %rms of 5 selected channels
rms_PC2_SbefD=cell2mat(resultsSbefD{7,2});
rms_PC3_SbefD=cell2mat(resultsSbefD{8,2});
rms_PC1_Dev=cell2mat(resultsDev{6,2}); 
rms_PC2_Dev=cell2mat(resultsDev{7,2});
rms_PC3_Dev=cell2mat(resultsDev{8,2});

filename_H=resultsSbefD{1,2};

for n=1:size(path,1)
    
currentpath=cell2mat(path(n,:));
cd(currentpath)

if size(filename_H,1)~=1
    currentfile=cell2mat(filename_H(n,:));
    filename_short=currentfile(1,1:18);
else
    currentfile=filename_H;
    filename_short=currentpath(1,68:87);
end


%check all different combinations of PC 1, 2 & 3
if PCA_check_SbefD(n)==0
        
        disp(['File not used: ' filename_short])
        continue
        
    elseif PCA_check_SbefD(n)==1 && PCA_check_Dev(n)==1
    
           standard=rms_PC1_SbefD(n,:); %rms over 5 selected channels related to PC1 for standard
           
           deviant=rms_PC1_Dev(n,:); %rms over 5 selected channels related to PC1 for Deviant
   
    elseif PCA_check_SbefD(n)==1 && PCA_check_Dev(n)==2

           standard=rms_PC1_SbefD(n,:); %rms over 5 selected channels related to PC1 for standard
 
           deviant=rms_PC2_Dev(n,:); %rms over 5 selected channels related to PC2 for Deviant
           
     elseif PCA_check_SbefD(n)==1 && PCA_check_Dev(n)==3

           standard=rms_PC1_SbefD(n,:); %rms over 5 selected channels related to PC1 for standard

           deviant=rms_PC3_Dev(n,:); %rms over 5 selected channels related to PC3 for Deviant
           
    elseif PCA_check_SbefD(n)==2 && PCA_check_Dev(n)==1

           standard=rms_PC2_SbefD(n,:); %rms over 5 selected channels related to PC2 for standard

           deviant=rms_PC1_Dev(n,:); %rms over 5 selected channels related to PC1 for Deviant  
           
    elseif PCA_check_SbefD(n)==2 && PCA_check_Dev(n)==2

           standard=rms_PC2_SbefD(n,:); %rms over 5 selected channels related to PC2 for standard

           deviant=rms_PC2_Dev(n,:); %rms over 5 selected channels related to PC2 for Deviant
           
     elseif PCA_check_SbefD(n)==2 && PCA_check_Dev(n)==3

           standard=rms_PC2_SbefD(n,:); %rms over 5 selected channels related to PC2 for standard

           deviant=rms_PC3_Dev(n,:); %rms over 5 selected channels related to PC3 for Deviant
           
     elseif PCA_check_SbefD(n)==3 && PCA_check_Dev(n)==1

           standard=rms_PC3_SbefD(n,:); %rms over 5 selected channels related to PC3 for standard

           deviant=rms_PC1_Dev(n,:); %rms over 5 selected channels related to PC1 for Deviant
           
     elseif PCA_check_SbefD(n)==3 && PCA_check_Dev(n)==2

           standard=rms_PC3_SbefD(n,:); %rms over 5 selected channels related to PC3 for standard

           deviant=rms_PC2_Dev(n,:); %rms over 5 selected channels related to PC2 for Deviant   
           
     elseif PCA_check_SbefD(n)==3 && PCA_check_Dev(n)==3

           standard=rms_PC3_SbefD(n,:); %rms over 5 selected channels related to PC3 for standard

           deviant=rms_PC3_Dev(n,:); %rms over 5 selected channels related to PC3 for Deviant
           
end  
           
mmseq=deviant-standard;


tt4=1800+tbeft; %tone onset + baseline
teval=1200; % time after tone that should be evaluated

sl=1/sf*1000;
tp=round(tt4/sl); %timepoint in samples 

tlength=round(teval/sl);
tone=mmseq(:,tp:tp+tlength); % time interval that shall be evaluated
rms_tone=rms(tone,1);

[pks, ptm]=findpeaks(rms_tone, 'MinPeakHeight', 5e-15); %search for mismatch peaks
ptm_t=ptm.*sl; %timing of peaks
[col1]=find(ptm_t>= 100); %only consider peaks after 100ms
peak_pos=[ptm_t',pks'];
sel_peaks=peak_pos(col1,:);
mm_peaks{n,:}=sel_peaks; %save all peaks in cell

dev_tone=deviant(:,tp:tp+tlength); %create tones individually to save sequences
std_tone=standard(:,tp:tp+tlength);

dev{n,:}=dev_tone;
stand{n,:}=std_tone;

dev_whole{n,:}=deviant;
std_whole{n,:}=standard;


%plot
    for j= 1:size(tone,2)
        time(j)=sl*j;
    end

figure
plot(time, rms_tone, 'LineWidth', 2)
hold on
plot(time, tone, 'LineWidth', 2)
hold on
plot(sel_peaks(:,1), sel_peaks(:,2), 'ok', 'MarkerFaceColor', 'k')
xlim([0 teval])
xlabel('Time in ms','FontSize',18)
ylabel('Magnetic field strength in Tesla','FontSize',18)
set(gca,'fontsize',14)
box off
saveas(figure(1), [filename_short, '_MM_between'], 'bmp')
close all
end

index=find(PCA_check_SbefD(:,1)~=0);

if size(filename_H,1)~=1
    name=filename_H(index,:);
else
    name=filename_H;
end

relevant_mm_peaks=mm_peaks(index,:);
dev_t=dev(index,:);
std_t=stand(index,:);

dev_w=dev_whole(index,:);
std_w=std_whole(index,:);

MMbetween={'filename', name; 'MM_peaks', relevant_mm_peaks;...
    'data_standard_whole', std_w; 'data_deviant_whole', dev_w;...
    'data_standard_t4', std_t; 'data_deviant_t4', dev_t};

end