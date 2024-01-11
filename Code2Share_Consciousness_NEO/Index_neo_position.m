function [pos_indS, pos_indD]=Index_neo_position(Top5)
%1. calculate mean position for each principal component
%fetal position ind of 1 means BEL or QL and 0 meand SL

for ki=1:size(Top5{1,2},2)
    input=cell2mat(Top5{1,2}{1,ki}(:,2:end));
    mean_pos=mean(input,1);
    mean_chan_pos_sPC1(ki,:)=mean_pos;
end

for ki=1:size(Top5{2,2},2)
    input=cell2mat(Top5{2,2}{1,ki}(:,2:end));
    mean_pos=mean(input,1);
    mean_chan_pos_sPC2(ki,:)=mean_pos;
end

for ki=1:size(Top5{3,2},2)
    input=cell2mat(Top5{3,2}{1,ki}(:,2:end));
    mean_pos=mean(input,1);
    mean_chan_pos_sPC3(ki,:)=mean_pos;
end

for ki=1:size(Top5{4,2},2)
    input=cell2mat(Top5{4,2}{1,ki}(:,2:end));
    mean_pos=mean(input,1);
    mean_chan_pos_dPC1(ki,:)=mean_pos;
end

for ki=1:size(Top5{5,2},2)
    input=cell2mat(Top5{5,2}{1,ki}(:,2:end));
    mean_pos=mean(input,1);
    mean_chan_pos_dPC2(ki,:)=mean_pos;
end

for ki=1:size(Top5{6,2},2)
    input=cell2mat(Top5{6,2}{1,ki}(:,2:end));
    mean_pos=mean(input,1);
    mean_chan_pos_dPC3(ki,:)=mean_pos;
end
%2. comparison with neonatal position
for io=1:size(mean_chan_pos_sPC1,1)
        if mean_chan_pos_sPC1(io,2)<-15 
        index=0;
        elseif mean_chan_pos_sPC1(io,2)>15
        index=0;
        elseif mean_chan_pos_sPC1(io,1)>15
        index=0;
        elseif mean_chan_pos_sPC1(io,1)<-15
        index=0;
        else
        index=1;
        end
        pos_ind(1,io)=index;
end

for io=1:size(mean_chan_pos_sPC2,1)
        if mean_chan_pos_sPC2(io,2)<-15 
        index=0;
        elseif mean_chan_pos_sPC2(io,2)>15
        index=0;
         elseif mean_chan_pos_sPC2(io,1)>15
        index=0;
        elseif mean_chan_pos_sPC2(io,1)<-15
        index=0;
        else
        index=1;
        end
        pos_ind(2,io)=index;
end

for io=1:size(mean_chan_pos_sPC3,1)
        if mean_chan_pos_sPC3(io,2)<-15 
        index=0;
        elseif mean_chan_pos_sPC3(io,2)>15
        index=0;
        elseif mean_chan_pos_sPC3(io,1)>15
        index=0;
        elseif mean_chan_pos_sPC3(io,1)<-15
        index=0;
        else
        index=1;
        end
        pos_ind(3,io)=index;
end

for io=1:size(mean_chan_pos_dPC1,1)
        if mean_chan_pos_dPC1(io,2)<-15 
        index=0;
        elseif mean_chan_pos_dPC1(io,2)>15
        index=0;
        elseif mean_chan_pos_dPC1(io,1)>15
        index=0;
        elseif mean_chan_pos_dPC1(io,1)<-15
        index=0;
        else
        index=1;
        end
        pos_ind(4,io)=index;
end

for io=1:size(mean_chan_pos_dPC2,1)
        if mean_chan_pos_dPC2(io,2)<-15 
        index=0;
        elseif mean_chan_pos_dPC2(io,2)>15
        index=0;
        elseif mean_chan_pos_dPC2(io,1)>15
        index=0;
        elseif mean_chan_pos_dPC2(io,1)<-15
        index=0;
        else
        index=1;
        end
        pos_ind(5,io)=index;
end
for io=1:size(mean_chan_pos_dPC3,1)
        if mean_chan_pos_dPC3(io,2)<-15 
        index=0;
        elseif mean_chan_pos_dPC3(io,2)>15
        index=0;
         elseif mean_chan_pos_dPC3(io,1)>15
        index=0;
        elseif mean_chan_pos_dPC3(io,1)<-15
        index=0;
        else
        index=1;
        end
        pos_ind(6,io)=index;
end
pos_indS=pos_ind(1:3,:);
pos_indD=pos_ind(4:6,:);
end