function []=JM_peak_location_plot_v2(path, sf, trig, tbeft, taft, results, PCA_check)
%this function creates plots to display rms trace and peaks for each subject
%in comparison to first version, this function only marks the peaks, that
%fulfill plusminus noise criterion.
 % function for finding selected peaks within the vector that the peak
% search created. Peak time is transformend from samples to seconds and
% peaks are searched in the time window 100-400ms after each tone onset and
% 100-900ms after the last tone. 

trial_data_all=results{5,2}; %trial data from all channels
trial_data_high5=results{6,2}; %5 selected channels
trial_data_high52=results{7,2};
trial_data_high53=results{8,2};
sensors_high5=results{9,2}; %5 selected channels
sensors_high52=results{10,2};
sensors_high53=results{11,2};
if size(results,1)>32
selected_peaks_plusminus=results{31,2};
selected_peaks_plusminus2=results{32,2};
selected_peaks_plusminus3=results{33,2};
else
selected_peaks_plusminus=results{27,2};
selected_peaks_plusminus2=results{28,2};
selected_peaks_plusminus3=results{29,2};    
end
filename_H=results{1,2};


 sl=1/sf*1000;
 for j= 1:length(trial_data_all{1})
    time(j)=sl*j-tbeft;
 end


for n=1:size(PCA_check,1) 

    %set path and filename
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
 % plots for all peaks including traces of 5 highest sensors and rms PC1

    
    current_peaks=selected_peaks_plusminus{n}; %take peaks with plusminus noise criterion
    current_sensors=sensors_high5{n};
    current_rms=trial_data_high5{n};
    current_data=trial_data_all{n};
    figure
    plot([0,0],[min(min(current_data(current_sensors, :))), max(max(current_data(current_sensors, :)))+1e-15], 'b', 'LineWidth',2)
    hold on
    plot([600,600],[min(min(current_data(current_sensors, :))), max(max(current_data(current_sensors, :)))+1e-15], 'b', 'LineWidth',2)
    hold on
    plot([1200,1200],[min(min(current_data(current_sensors, :))), max(max(current_data(current_sensors, :)))+1e-15], 'b', 'LineWidth',2)
    hold on
    plot([1800,1800],[min(min(current_data(current_sensors, :))), max(max(current_data(current_sensors, :)))+1e-15], 'b', 'LineWidth',2)
    hold on
    plot(time, current_rms, 'r', 'LineWidth', 2.5)
    hold on
    plot(time, current_data(current_sensors, :), 'k')
    hold on
    plot(current_peaks(:,1), current_peaks(:,2), 'ok', 'MarkerFaceColor', 'k')
    %title('Peaks in selected time windows - Orig, PC1 - DevD')
    xlim([-tbeft taft])
    xlabel('Time in ms','FontSize',18)
    ylabel('Magnetic field strength in Tesla','FontSize',18)
    set(gca,'fontsize',14)
    box off
    set(gcf, 'Position', [1050, 650, 1050, 650])
    saveas(figure(1), [filename_short, '_Orig_PC1_', trig], 'bmp')
    close all


elseif PCA_check(n)==2
 % plots for all peaks including traces of 5 highest sensors and rms PC2

    current_peaks=selected_peaks_plusminus2{n};
    current_sensors=sensors_high52{n};
    current_rms=trial_data_high52{n};
    current_data=trial_data_all{n};
    figure
    plot([0,0],[min(min(current_data(current_sensors, :))), max(max(current_data(current_sensors, :)))+1e-15], 'b', 'LineWidth',2)
    hold on
    plot([600,600],[min(min(current_data(current_sensors, :))), max(max(current_data(current_sensors, :)))+1e-15], 'b', 'LineWidth',2)
    hold on
    plot([1200,1200],[min(min(current_data(current_sensors, :))), max(max(current_data(current_sensors, :)))+1e-15], 'b', 'LineWidth',2)
    hold on
    plot([1800,1800],[min(min(current_data(current_sensors, :))), max(max(current_data(current_sensors, :)))+1e-15], 'b', 'LineWidth',2)
    hold on
    plot(time, current_rms, 'r', 'LineWidth', 2.5)
    hold on
    plot(time, current_data(current_sensors, :), 'k')
    hold on
    plot(current_peaks(:,1), current_peaks(:,2), 'ok', 'MarkerFaceColor', 'k')
    %title('Peaks in selected time windows - Orig, PC2 - DevD')
    xlim([-tbeft taft])
    xlabel('Time in ms','FontSize',18)
    ylabel('Magnetic field strength in Tesla','FontSize',18)
    set(gca,'fontsize',14)
    box off
    set(gcf, 'Position', [1050, 650, 1050, 650])
    saveas(figure(1), [filename_short, '_Orig_PC2_', trig], 'bmp')
    close all
    
    elseif PCA_check(n)==3
 % plots for all peaks including traces of 5 highest sensors and rms PC3

    current_peaks=selected_peaks_plusminus3{n};
    current_sensors=sensors_high53{n};
    current_rms=trial_data_high53{n};
    current_data=trial_data_all{n};
    figure
    plot([0,0],[min(min(current_data(current_sensors, :))), max(max(current_data(current_sensors, :)))+1e-15], 'b', 'LineWidth',2)
    hold on
    plot([600,600],[min(min(current_data(current_sensors, :))), max(max(current_data(current_sensors, :)))+1e-15], 'b', 'LineWidth',2)
    hold on
    plot([1200,1200],[min(min(current_data(current_sensors, :))), max(max(current_data(current_sensors, :)))+1e-15], 'b', 'LineWidth',2)
    hold on
    plot([1800,1800],[min(min(current_data(current_sensors, :))), max(max(current_data(current_sensors, :)))+1e-15], 'b', 'LineWidth',2)
    hold on
    plot(time, current_rms, 'r', 'LineWidth', 2.5)
    hold on
    plot(time, current_data(current_sensors, :), 'k')
    hold on
    plot(current_peaks(:,1), current_peaks(:,2), 'ok', 'MarkerFaceColor', 'k')
    %title('Peaks in selected time windows - Orig, PC2 - DevD')
    xlim([-tbeft taft])
    xlabel('Time in ms','FontSize',18)
    ylabel('Magnetic field strength in Tesla','FontSize',18)
    set(gca,'fontsize',14)
    box off
    set(gcf, 'Position', [1050, 650, 1050, 650])
    saveas(figure(1), [filename_short, '_Orig_PC3_', trig], 'bmp')
    close all
end

end
end