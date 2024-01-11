function [MMwithin]=MMwithinSequence_v3(path, sf, tbeft, trig, results, PCA_check)
%v3 uses rms of 5 channels instead of mean


rms_PC1=cell2mat(results{6,2}); %rms 5 selected channels
rms_PC2=cell2mat(results{7,2});
rms_PC3=cell2mat(results{8,2});
filename_H=results{1,2};

%use PCA_check to distinguish between files that are not used (PCA_check=0)
%where PC1 is used (PCA_check=1) and where PC2 is used (PCA_check=2)

for n=1:size(PCA_check,1)
    
currentpath=cell2mat(path(n,:));
cd(currentpath)
if size(filename_H,1)~=1
    currentfile=cell2mat(filename_H(n,:));
    filename_short=currentfile(1,1:18);
else
    currentfile=filename_H;
    filename_short=currentpath(1,68:87);
end
    

    if PCA_check(n)==0
        
        disp(['File not used: ' filename_short])
        continue
        
    elseif PCA_check(n)==1

        current_data=rms_PC1(n,:); %rms over 5 selected channels related to PC1


    elseif PCA_check(n)==2 

        current_data=rms_PC2(n,:); %rms over 5 selected channels related to PC2
    
    elseif PCA_check(n)==3 

        current_data=rms_PC3(n,:); %rms over 5 selected channels related to PC2
    end
    
tt3=1200+tbeft; %tone onset + baseline
tt4=1800+tbeft; %tone onset + baseline
teval=600; % time after each tone that should be evaluated

sl=1/sf*1000;
t1=round(tt3/sl); %timepoint in samples 
t2=round(tt4/sl);
tlength=round(teval/sl);
tone3=current_data(:,t1:t1+tlength); % time interval that shall be substracted
tone4=current_data(:,t2:t2+tlength);
mmseq=tone4-tone3;
rms_mmseq=rms(mmseq,1);
[pks, ptm]=findpeaks(rms_mmseq, 'MinPeakHeight', 5e-15); %search for mismatch peaks
ptm_t=ptm.*sl; %timing of peaks
[col1]=find(ptm_t>= 100); %only consider peaks after 100ms
peak_pos=[ptm_t',pks'];
sel_peaks=peak_pos(col1,:);
mm_peaks{n,:}=sel_peaks; %save all peaks in cell
data_tone3{n,:}=tone3;
data_tone4{n,:}=tone4;

%plot
    for j= 1:size(mmseq,2)
        time(j)=sl*j;
    end

figure
plot(time, rms_mmseq, 'LineWidth', 2)
hold on
plot(time, mmseq,  'LineWidth', 2)
hold on
plot(sel_peaks(:,1), sel_peaks(:,2), 'ok', 'MarkerFaceColor', 'k')
xlim([0 teval])
xlabel('Time in ms','FontSize',18)
ylabel('Magnetic field strength in Tesla','FontSize',18)
set(gca,'fontsize',14)
box off
saveas(figure(1), [filename_short, '_MM_Seq_', trig], 'bmp')
close all



end

index=find(PCA_check(:,1)~=0);
if size(filename_H,1)~=1
    name=filename_H(index,:);
else
    name=filename_H;
end

relevant_mm_peaks=mm_peaks(index,:);
d_tone3=data_tone3(index,:);
d_tone4=data_tone4(index,:);

MMwithin={'filename', name; 'MM_peaks', relevant_mm_peaks; 'data_tone3', d_tone3; 'data_tone4', d_tone4};

end