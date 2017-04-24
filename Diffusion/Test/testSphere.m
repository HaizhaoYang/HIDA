%    Points = set of N points on the sphere
%    Tree   = diffusion wavelet packet tree for the diffusion operator
%    F      = first compression example function
%    G      = second compression example function

clc;
close all;
clear;           % we'll have more than enough data as it is
%% set up parameters
N = 2000        % number of points to use
MBest = 400         % best M term approximation
M = 400        % precomputed M term appproximation
nvars = M;%N-1;%11  % number of leading eigenvalues used in the diffusion map
numNearest = 10; % number of nearest points


%% construct diffusion operators and representations
if MBest<=M
    filename = sprintf('./Data/sphere_N_%d_M_%d.mat', N,M);
else
    filename = sprintf('./Data/sphere_N_%d_M_%d.mat', N,MBest);
end

% check to see if the data has already been generated
if exist(filename)
    fprintf('loading %s data file ...', filename);
    load(filename);
    fprintf('\n');
else
    
    fprintf('could not find %s, generating data from scratch ... \n', filename);
    Points = GeneratePoints('Sphere', N);
    
    % build the two functions used in the compression examples
    [Theta,Phi,Rho] = cart2sph(Points(:,1), Points(:,2), Points(:,3));
    
    % F is trigonemtric poly + exponential at different frequncies
    %   F = sin(3*Theta.^8) .* exp(-5*Phi.^2);
    % build F
    F = sin(Theta);	% background
    
    % add swirl
    Pt = Points(26,:);
    
    atria = nn_prepare(Points);
    [count,neighbors] = range_search(Points, atria, Pt, 1.0);
    
    Delta = 10;
    idxs = neighbors{1};
    dists = neighbors{2};
    
    ips = Points(idxs, :)*Pt';
    
    %F1(neighbors{1}(j)) = Points(neighbors{1}(j),1)^4*sin(20*ip) * exp(-Delta*neighbors{2}(j)^2);
    F(idxs) = F(idxs)+Points(idxs,1).^4.*sin(20*ips).*exp(-Delta*dists.^2)';
    
    idxs = find(0 < Theta & Theta < pi/4 & 0 < Phi & Phi < pi/6);
    F(idxs) = -1.0;
    
    G = zeros(N,1);
    G = G+.5*cos(.6*Theta).*exp(-3*Phi);
    
    %G = cos(Phi).^3-sin(Theta.^2)+cos(3*Phi-Theta);
    %idxs = find(-pi/2 <= Phi & Phi <= pi/2 & pi/4<= Theta & Theta <= pi/2);
    %G(idxs) = 1.4;
    
    % construct a diffusion operator and an oversampled diffusion operator
    T = MakeDiffusion(Points, numNearest, struct('Delta', 160, 'Threshold', 1e-3, 'Normalization', 'beltrami'));
    
    % compute eigenvalues and vectors for the oversampled sphere
    fprintf('solving eigenproblem for T ... '); TIME = cputime;
    [V, D] = eigs(T,nvars);
    % flip the eigenvectors
    V = fliplr(V);
    basisDM = V;
    dmTime = cputime - TIME;
    
    % VBDM
    fprintf('computing VBDM ... '); TIME = cputime;
    k = numNearest;%512;
    k2 = numNearest;%64;
    operator = 2;
    
    [basisWeightedDM,eigs3,epsilon3,~,peq3,qest3,dim3] = VBDM(Points,k,k2,nvars,operator);
    wdmTime = cputime - TIME;
    
    % now compute the tree
    fprintf('constructing diffusion wavelets ... '); TIME = cputime;
    Tree = DWPTree(T, 12, 1e-5, struct('StopFcns', 10, 'ExtendBases', false)); % Haizhao: check why 14 levels doesn't work!
    dwTime = cputime - TIME;
    
    fprintf('\n\nsaving sphere data to %s ...', filename);
    SphereDataLoaded = 1;
    save(filename,'-v7.3');
end


fprintf('\n');
fprintf('Time for diffusion map                   %g seconds\n', dmTime);
fprintf('Time for weighted diffusion map          %g seconds\n', wdmTime);
fprintf('Time for diffusion wavelet               %g seconds\n', dwTime);

%% test 1
% compute the coefficients by the diffusion map
coeffsDM = basisDM'*G; % I think Mauro's example is wrong and we should use basisDM'*G
[junk idxs] = sort(abs(coeffsDM),'descend');
coeffsDM(idxs(MBest+1:nvars)) = 0.0;
ReconG1 = basisDM*coeffsDM;

% compute the coefficients by the weighted diffusion map
ww = peq3./qest3/N;
coeffsWeightedDM = basisWeightedDM'*(ww.*G);
[junk,idxs] = sort(abs(coeffsWeightedDM),'descend');
coeffsDM(idxs(MBest+1:nvars)) = 0.0;
ReconG2 = basisWeightedDM*coeffsWeightedDM;

% compute the coefficients by the diffusion wavelet
coeffsDW = DWCoeffs(Tree, G);
Best   = DWBest(coeffsDW);
BestList = DWUnpack(Best);
if 1
    [junk,idxs] = sort(abs(BestList(:,4)),'descend');
    BestList(idxs(MBest+1:end),4) = 0.0;
    BestM = BestList;
else
    BestM = BestList(1:MBest,:);
end
Best = DWPack(Tree, BestM);
ReconG3 = DWRecon(Tree, Best);

figure;
subplot(2,2,1);DrawSphereFcn(Points, G);title('Original figure');
subplot(2,2,2);DrawSphereFcn(Points, ReconG1);
head = sprintf('Recon. by DM, rel. err.: %f',norm(G-ReconG1)/norm(G)); title(head);
subplot(2,2,3);DrawSphereFcn(Points, ReconG2);
head = sprintf('Recon. by wDM, rel. err.: %f',norm(G-ReconG2)/norm(G)); title(head);
subplot(2,2,4);DrawSphereFcn(Points, full(ReconG3));
head = sprintf('Recon. by DW, rel. err.: %f',norm(G-ReconG3)/norm(G)); title(head);

fprintf('\n');
fprintf('Figure1: Reconstruction of G from the top %d eigenfunction coefficients.\n',M);
fprintf('Figure2: Reconstruction of G from the top %d weighted eigenfunction coefficients.\n',M);
fprintf('Figure3: Reconstruction of G from the top %d diffusion wavelet coefficients.\n',M);

fprintf('\n');
fprintf('Reconstruction accuracy of diffusion map               = %g\n', norm(G-ReconG1)/norm(G));
fprintf('Reconstruction accuracy of weighted diffusion map      = %g\n', norm(G-ReconG2)/norm(G));
fprintf('Reconstruction accuracy of diffusion wavelet           = %g\n', norm(G-ReconG3)/norm(G));

%% test 2

% compute the coefficients by the diffusion map
coeffsDM = basisDM'*F; % I think Mauro's example is wrong and we should use basisDM'*G
[junk idxs] = sort(abs(coeffsDM),'descend');
coeffsDM(idxs(MBest+1:nvars)) = 0.0;
ReconF1 = basisDM*coeffsDM;

% compute the coefficients by the weighted diffusion map
ww = peq3./qest3/N;
coeffsWeightedDM = basisWeightedDM'*(ww.*F);
[junk,idxs] = sort(abs(coeffsWeightedDM),'descend');
coeffsDM(idxs(MBest+1:nvars)) = 0.0;
ReconF2 = basisWeightedDM*coeffsWeightedDM;

% compute the coefficients by the diffusion wavelet
coeffsDW = DWCoeffs(Tree, F);
Best   = DWBest(coeffsDW);
BestList = DWUnpack(Best);
if 1
    [junk,idxs] = sort(abs(BestList(:,4)),'descend');
    BestList(idxs(MBest+1:end),4) = 0.0;
    BestM = BestList;
else
    BestM = BestList(1:MBest,:);
end
Best = DWPack(Tree, BestM);
ReconF3 = DWRecon(Tree, Best);

figure;
subplot(2,2,1);DrawSphereFcn(Points, F);title('Original figure');
subplot(2,2,2);DrawSphereFcn(Points, ReconF1);
head = sprintf('Recon. by DM, rel. err.: %f',norm(F-ReconF1)/norm(F)); title(head);
subplot(2,2,3);DrawSphereFcn(Points, ReconF2);
head = sprintf('Recon. by wDM, rel. err.: %f',norm(F-ReconF2)/norm(F)); title(head);
subplot(2,2,4);DrawSphereFcn(Points, full(ReconF3));
head = sprintf('Recon. by DW, rel. err.: %f',norm(F-ReconF3)/norm(F)); title(head);


fprintf('\n');
fprintf('Figure1: Reconstruction of F from the top %d eigenfunction coefficients.\n',M);
fprintf('Figure2: Reconstruction of F from the top %d weighted eigenfunction coefficients.\n',M);
fprintf('Figure3: Reconstruction of F from the top %d diffusion wavelet coefficients.\n',M);

fprintf('\n');
fprintf('Reconstruction accuracy of diffusion map               = %g\n', norm(F-ReconF1)/norm(F));
fprintf('Reconstruction accuracy of weighted diffusion map      = %g\n', norm(F-ReconF2)/norm(F));
fprintf('Reconstruction accuracy of diffusion wavelet           = %g\n', norm(F-ReconF3)/norm(F));

