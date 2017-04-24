% Figure14.m
%
% Generate Figure 14 in section 7.3.
%
% Original code by Bremer et al.
%
% Modified by Haizhao Yang to replace some old routines, results may 
% differ from those in the paper by Bremer et al.

if ~exist('SphereDataLoaded')
   fprintf('Figure14.m: Sphere data not found; run "SphereData.m" first\n');
   return;
end

M = 200;

Coeffs = DWCoeffs(Tree, G);
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

ReconG1 = DWRecon(Tree, Best);


EigCoeffs = Eigs*G;
[junk idxs] = sort(abs(EigCoeffs),'descend');
EigCoeffs(idxs(M+1:N)) = 0.0;
ReconG2 = Eigs'*EigCoeffs;


% compute the coefficients by the diffusion map of the small data set
[V, D] = eigs(T,N-1);
V = fliplr(V);
basisDM = V;
coeffsDM = basisDM'*G; % I think Mauro's example is wrong and we should use basisDM'*G
[junk idxs] = sort(abs(coeffsDM),'descend');
coeffsDM(idxs(M+1:end)) = 0.0;
ReconG3 = basisDM*coeffsDM;

pic = figure;set(pic, 'Position', [200, 200, 1200, 1200]);
subplot(2,2,1);DrawSphereFcn(Points, G);title('Original figure');
subplot(2,2,2);DrawSphereFcn(Points, full(ReconG1));
head = sprintf('Recon. by DW, rel. err.: %f',norm(G-ReconG1)/norm(G)); title(head);
subplot(2,2,3);DrawSphereFcn(Points, ReconG2);
head = sprintf('Recon. by DM by Mauro, rel. err.: %f',norm(G-ReconG2)/norm(G)); title(head);
subplot(2,2,4);DrawSphereFcn(Points, ReconG3);
head = sprintf('Recon. by DM by ours, rel. err.: %f',norm(G-ReconG3)/norm(G)); title(head);
axesHandles = findobj(get(pic,'Children'), 'flat','Type','axes');
axis(axesHandles,'square')