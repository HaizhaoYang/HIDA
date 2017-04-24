% Figure11.m
% Original code by Bremer et al.
%
% Generate Figure 11 in section 7.2.
%
% Modified by Haizhao Yang to replace some old routines, results may 
% differ from those in the paper by Bremer et al.

if ~exist('SphereDataLoaded')
    fprintf('Figure6.m: Sphere data not found; run "SphereData.m" first\n');
    return;
end

M = 200;

close all;

Coeffs = DWCoeffs(Tree, F);
Best   = DWBest(Coeffs);
BestList = DWUnpack(Best);
if 1 % Method by Haizhao Yang
    [junk,idxs] = sort(abs(BestList(:,4)),'descend');
    BestList(idxs(M+1:end),4) = 0.0;
    BestM = BestList;
else % Method by Bremer et al.
    BestM = BestList(1:M,:);
end
Best = DWPack(Tree, BestM);

ReconF1 = DWRecon(Tree, Best);


EigCoeffs = Eigs*F;
[junk idxs] = sort(abs(EigCoeffs),'descend');
EigCoeffs(idxs(M+1:N)) = 0.0;
ReconF2 = Eigs'*EigCoeffs;


% compute the coefficients by the diffusion map of the small data set
[V, D] = eigs(T,N-1);
V = fliplr(V);
basisDM = V;
coeffsDM = basisDM'*F; % I think Mauro's example is wrong and we should use basisDM'*G
[junk idxs] = sort(abs(coeffsDM),'descend');
coeffsDM(idxs(M+1:end)) = 0.0;
ReconF3 = basisDM*coeffsDM;

pic = figure;set(pic, 'Position', [200, 200, 1200, 1200]);
subplot(2,2,1);DrawSphereFcn(Points, F);title('Original figure');
subplot(2,2,2);DrawSphereFcn(Points, full(ReconF1));
head = sprintf('Recon. by DW, rel. err.: %f',norm(F-ReconF1)/norm(F)); title(head);
subplot(2,2,3);DrawSphereFcn(Points, ReconF2);
head = sprintf('Recon. by DM by Mauro, rel. err.: %f',norm(F-ReconF2)/norm(F)); title(head);
subplot(2,2,4);DrawSphereFcn(Points, ReconF3);
head = sprintf('Recon. by DM by ours, rel. err.: %f',norm(F-ReconF3)/norm(F)); title(head);
axesHandles = findobj(get(pic,'Children'), 'flat','Type','axes');
axis(axesHandles,'square')