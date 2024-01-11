function [PCA_check]=PCA_check_creation_function(sel_comp, Neo_pos_index)
sel_comp=sel_comp';
%create count of included datasets
dat_count=[1:size(Neo_pos_index,2)]';
%determine positions that need to be filled with zeros
[A,b]=ismember(dat_count, sel_comp(:,1));
b(A)=sel_comp(:,2);
PCA_check=b;
end