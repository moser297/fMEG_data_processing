function [idxS, idxD, within_dist]=Index_distance_within_chan_clusters_neo(Top5)

for ki=1:size(Top5{1,2},2)
    input=cell2mat(Top5{1,2}{1,ki}(:,2:end));
    within_dist=squareform(pdist(input));
    mask = tril(true(size(within_dist)),-1);
    within_d_m=mean(within_dist(mask));
    %within_d_max=max(within_dist(mask));
    mean_within_dist_sPC1(ki,:)=within_d_m;
    %max_within_dist_PC1(ki,:)=within_d_max;
end

for ki=1:size(Top5{2,2},2)
    input=cell2mat(Top5{2,2}{1,ki}(:,2:end));
    within_dist=squareform(pdist(input));
    mask = tril(true(size(within_dist)),-1);
    within_d_m=mean(within_dist(mask));
    %within_d_max=max(within_dist(mask));
    mean_within_dist_sPC2(ki,:)=within_d_m;
    %max_within_dist_PC2(ki,:)=within_d_max;
end

for ki=1:size(Top5{3,2},2)
    input=cell2mat(Top5{3,2}{1,ki}(:,2:end));
    within_dist=squareform(pdist(input));
    mask = tril(true(size(within_dist)),-1);
    within_d_m=mean(within_dist(mask));
    %within_d_max=max(within_dist(mask));
    mean_within_dist_sPC3(ki,:)=within_d_m;
    %max_within_dist_PC3(ki,:)=within_d_max;
end

for ki=1:size(Top5{4,2},2)
    input=cell2mat(Top5{4,2}{1,ki}(:,2:end));
    within_dist=squareform(pdist(input));
    mask = tril(true(size(within_dist)),-1);
    within_d_m=mean(within_dist(mask));
    %within_d_max=max(within_dist(mask));
    mean_within_dist_dPC1(ki,:)=within_d_m;
    %max_within_dist_PC1(ki,:)=within_d_max;
end

for ki=1:size(Top5{5,2},2)
    input=cell2mat(Top5{5,2}{1,ki}(:,2:end));
    within_dist=squareform(pdist(input));
    mask = tril(true(size(within_dist)),-1);
    within_d_m=mean(within_dist(mask));
    %within_d_max=max(within_dist(mask));
    mean_within_dist_dPC2(ki,:)=within_d_m;
    %max_within_dist_PC1(ki,:)=within_d_max;
end

for ki=1:size(Top5{6,2},2)
    input=cell2mat(Top5{6,2}{1,ki}(:,2:end));
    within_dist=squareform(pdist(input));
    mask = tril(true(size(within_dist)),-1);
    within_d_m=mean(within_dist(mask));
    %within_d_max=max(within_dist(mask));
    mean_within_dist_dPC3(ki,:)=within_d_m;
    %max_within_dist_PC1(ki,:)=within_d_max;
end

within_dist={ 'sPC1', mean_within_dist_sPC1; ...
    'sPC2', mean_within_dist_sPC2; ...
    'sPC3', mean_within_dist_sPC3;...
    'dPC1', mean_within_dist_dPC1;...
    'dPC2', mean_within_dist_dPC2;...
    'dPC3', mean_within_dist_dPC3};

for io=1:size(within_dist{1,2},1)
        if within_dist{1,2}(io,1)>10
        index=0;
        else
        index=1;
        end
        idxS(1,io)=index;
        if within_dist{2,2}(io,1)>10
        index=0;
        else
        index=1;
        end
        idxS(2,io)=index;
        if within_dist{3,2}(io,1)>10
        index=0;
        else
        index=1;
        end
        idxS(3,io)=index;
end
for io=1:size(within_dist{4,2},1)
        if within_dist{4,2}(io,1)>10
        index=0;
        else
        index=1;
        end
        idxD(1,io)=index;
        if within_dist{5,2}(io,1)>10
        index=0;
        else
        index=1;
        end
        idxD(2,io)=index;
        if within_dist{6,2}(io,1)>10
        index=0;
        else
        index=1;
        end
        idxD(3,io)=index;
end
end