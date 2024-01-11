function [PCA_check_S, PCA_check_D]=JM_check_matching_PCA_check(PCA_check_S, PCA_check_D)
%function checks, if PCA check was possitive for standard and deviant of a
%certain file. If it was not for one, the other is also set to 0
test=any(PCA_check_S,2); 
test2=any(PCA_check_D,2); 
mm=any(test-test2, 2);
PCA_check_S(mm)=0;
PCA_check_D(mm)=0;
end