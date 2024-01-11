function [Top5]=Chan_pos_top5_neo(resultsS, resultsD)
%extract information from all selected top 5 channels
channel_info=resultsS{2,2};
for j=1:size(channel_info,2)
SensSPC1=channel_info{1,j}(resultsS{9,2}{j,1},:);
SensSPC2=channel_info{1,j}(resultsS{10,2}{j,1},:);
SensSPC3=channel_info{1,j}(resultsS{11,2}{j,1},:);
SensDPC1=channel_info{1,j}(resultsD{9,2}{j,1},:);
SensDPC2=channel_info{1,j}(resultsD{10,2}{j,1},:);
SensDPC3=channel_info{1,j}(resultsD{11,2}{j,1},:);
SensPC1_w=channel_info{1,j}(resultsS{12,2}{j,1},:);
SensPC2_w=channel_info{1,j}(resultsS{13,2}{j,1},:);
SensPC3_w=channel_info{1,j}(resultsS{14,2}{j,1},:);
SensH=channel_info{1,j}(resultsS{15,2}{j,1},:);

Top5SPC1{j}=SensSPC1;
Top5SPC2{j}=SensSPC2;
Top5SPC3{j}=SensSPC3;
Top5DPC1{j}=SensDPC1;
Top5DPC2{j}=SensDPC2;
Top5DPC3{j}=SensDPC3;
Top5PC1_w{j}=SensPC1_w;
Top5PC2_w{j}=SensPC2_w;
Top5PC3_w{j}=SensPC3_w;
Top5H{j}=SensH;

end
Top5={'sPC1',Top5SPC1; 'sPC2',Top5SPC2; 'sPC3',Top5SPC3;...
    'dPC1',Top5DPC1; 'dPC2',Top5DPC2; 'dPC3',Top5DPC3;...
    'PC1_w', Top5PC1_w;'PC2_w', Top5PC2_w; 'PC3_w', Top5PC3_w; 'H',Top5H; };
end