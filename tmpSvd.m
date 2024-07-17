function [COH_permAbove,spSVmag_permAbove,COH_perm,spSVmag_perm] = tmpSvd(j2,allPerm,allPermN,V,K,F,COH,spSVmag)
% for each voxel, randomly permute tapers, using the same
% permutation across frequencies
[uPerm,sPerm,~] = pagesvd(...
    permute(  reshape(  j2(:,allPerm(randi(allPermN,V,1),:)' + (0:K:V*K-1))  ,[F K V])  ,[3 2 1])...
    ,'econ','vector');
COH_perm          = sPerm.^2./sum(sPerm.^2,1);
COH_permAbove     = COH_perm > COH;
spSVmag_perm      = abs(uPerm);
spSVmag_permAbove = spSVmag_perm > spSVmag;
