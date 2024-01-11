function data_filtered = KS_butterworthfilter(data, frequency, lower, upper)
% .................................... %
% by Katrin Sippel fMEG Tuebingen 2017 %
% 같같같같같같같같같같같같같같같같같같 %

data_filtered = zeros(size(data));

disp(strcat('... butterworth ', num2str(lower),'-', num2str(upper),'Hz'));%

for c = 1:1:size(data,1) 
    signal = data(c,:);
    [B,A] = butter(4,[lower upper]/(frequency/2));
    %sig_filtered = filter(B,A,signal-signal(1)); %only forward filtering
    sig_filtered = filtfilt(B,A,signal-signal(1)); %forward and backward filtering
    data_filtered(c,:) =sig_filtered;
end

end
