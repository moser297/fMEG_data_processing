function [selected_peaks]=FindPeaks(sf, peak_time, peaks)
% function for finding selected peaks within the vector that the peak
% search created. Peak time is transformend from samples to seconds and
% peaks are searched in the time window 100-400ms after each tone onset and
% 100-800ms after the last tone. Outputs selected peaks with time as unit,
% not samples. Delay from trigger to actual sound output of 27ms is already
% included as all triggers were shifted by 27ms

sl=1/sf*1000;
%figure
for i=1:size(peaks,2)
%Conversion from samples into time units
ptm_t=peak_time{i}.*sl-200;
% ptm2_t=peak_time2{i}.*sl-200;
% ptm_pc_t=peak_time_pc{i}.*sl-200;
% ptm_pc2_t=peak_time2_pc{i}.*sl-200;

pks=peaks{i};
% pks2=peaks2{i};
% pks_pc=peaks_pc{i};
% pks_pc2=peaks2_pc{i};

% figure
% plot(ptm_t,pks, 'r*')
% hold on
% plot(ptm2_t,pks2, 'b*')
% hold on
% plot(ptm_pc_t,pks_pc, 'm*')
% hold on
% plot(ptm_pc2_t,pks_pc2, 'c*')
% hold on
% plot([0,0],[0, 1e-14], 'r')
% hold on
% plot([600,600],[0, 1e-14], 'r')
% hold on
% plot([1200,1200],[0, 1e-14], 'r')
% hold on
% plot([1800,1800],[0, 1e-14], 'r')

% plot(ptm_t, pks, '*')
% hold on

%see which peaks are in range of 100-400ms after tone onset and 100-800ms
%for last tone
[col1]=find(ptm_t>= 100 & ptm_t<= 400);
[col2]=find(ptm_t>= 700 & ptm_t<= 1000);
[col3]=find(ptm_t>= 1300 & ptm_t<= 1600);
[col4]=find(ptm_t>= 1900 & ptm_t<= 2700);
peak_pos=[ptm_t',pks'];
sel_peaks=peak_pos([col1, col2, col3, col4],:);
selected_peaks{i}=sel_peaks;

end

end