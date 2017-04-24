function [psi,lambda,epsilon,peqoversample,peq,qest,dim] = VBDM(x,k,k2,nvars,operator,epsilon,dim)

%%% Inputs
%%% x       - N-by-n data set with N data points in R^n
%%% k       - number of nearest neighbors to use
%%% k2      - number of nearest neighbors to use to determine the "epsilon"
%%%             parameter
%%% nvars   - number of eigenfunctions/eigenvalues to compute
%%% operator- 1 - Laplace-Beltrami operator, 2 - generator of grad system
%%% dim     - optionally prescribe the intrinsic dimension of the manifold lying inside R^n
%%% epsilon - optionally choose an arbitrary "global" epsilon

%%% Outputs
%%% psi       - Eigenfunctions of the generator/Laplacian
%%% lambda       - Eigenvalues
%%% epsilon - scale, derived from the k2 nearest neighbors
%%% peqoversample - Invariant measure divided by the sampling measure
%%% peq     - Invariant measure
%%% qest    - Sampling measure
%%% dim     - estimated dimension


%%% Theory requires c2 = 1/2 - 2*alpha + 2*dim*alpha + dim*beta/2 + beta < 0
%%% The resulting operator will have c1 = 2 - 2*alpha + dim*beta + 2*beta
%%% Thus beta = (c1/2 - 1 + alpha)/(dim/2+1), since we want beta<0,
%%% natural choices are beta=-1/2 or beta = -1/(dim/2+1)

% Copyright 2017
% By John Harlin at Penn State Univ

N = size(x,1); %number of points

[d,inds] = knnCPU(x,x,k);

%%% Build ad hoc bandwidth function by autotuning epsilon for each pt.
epss = 2.^(-30:.1:10);
rho0 = sqrt(mean(d(:,2:k2).^2,2));

%%% Pre-kernel used with ad hoc bandwidth only for estimating dimension
%%% and sampling density
dt = d.^2./(repmat(rho0,1,k).*rho0(inds));

%%% Tune epsilon on the pre-kernel
dpreGlobal=zeros(1,length(epss));
for i=1:length(epss)
    dpreGlobal(i) = sum(sum(exp(-dt./(2*epss(i)))))/(N*k);
end

[maxval,maxind] = max(diff(log(dpreGlobal))./diff(log(epss)));

if (nargin < 7)
    dim=2*maxval
end

%%% Use ad hoc bandwidth function, rho0, to estimate the density
dt = exp(-dt./(2*epss(maxind)))/((2*pi*epss(maxind))^(dim/2));
dt = sparse(reshape(double(inds'),N*k,1),repmat(1:N,k,1),reshape(double(dt'),N*k,1),N,N,N*k)';
dt = (dt+dt')/2;

% sampling density estimate for bandwidth function
qest = (sum(dt,2))./(N*rho0.^(dim));

clear dt;

if (operator == 1)
    %%% Laplace-Beltrami, c1 = 0
    beta = -1/2;
    alpha = -dim/4 + 1/2;
elseif (operator == 2)
    %%% Kolmogorov backward operator, c1 = 1
    beta = -1/2;
    alpha = -dim/4;
end

c1 = 2 - 2*alpha + dim*beta + 2*beta;
c2=.5-2*alpha+2*dim*alpha+dim*beta/2+beta;

d = d.^2;

%%% bandwidth function rho(x) from the sampling density estimate
rho = qest.^(beta);
rho = rho/mean(rho);

%% construct the exponent of K^S_epsilon
d = d./repmat((rho),1,k);  % divide row j by rho(j)
d = d./rho(inds);

%%% Tune epsilon for the final kernel
if (nargin<6)
    for i=1:length(epss)
        s(i) = sum(sum(exp(-d./(4*epss(i))),2))/(N*k);
    end
    [~,maxind] = max(diff(log(s))./diff(log(epss)));
    epsilon = epss(maxind)
end

%%% K^S_epsilon with final choice of epsilon
d = exp(-d./(4*epsilon));
d = sparse(reshape(double(inds'),N*k,1),repmat(1:N,k,1),reshape(double(d'),N*k,1),N,N,N*k)';
clear inds;
d = (d+d')/2;   %%% symmetrize since this is the symmetric formulation

%%% q^S_epsilon (this is the sampling density estimate q(x) obtained from VB kernel)
qest = full((sum(d,2)./(rho.^dim)));

Dinv1 = spdiags(qest.^(-alpha),0,N,N);

%%% K = K^S_{epsilon,alpha}
d = Dinv1*d*Dinv1; % the "right" normalization

%%% S^2 =P*D, where P = diag(rho), D = q^S_{epsilon,alpha}
Ssquare = full((rho.^2).*(sum(d,2)));

%%% S^{-1}
Sinv = spdiags(Ssquare.^(-1/2),0,N,N);

%%% epsilon*Lhat +eye(I) = Sinv*K*Sinv - P^{-2} + eye(I)
d = Sinv*d*Sinv - spdiags(rho.^(-2)-1,0,N,N); %%% "left" normalization

opts.maxiter = 200;

if (nvars > 0)
    
    [psi,lambda] = eigs(d,nvars,1,opts);
    
    % this is eigenvalue of lambda = eig(eye(I) + epsilon*Lhat)
    % Since lambda^{1/epsilon} ---> e^(eig(Lhat))
    % and eig(Lhat) = log(lambda^{1/epsilon})
    lambda = (diag(lambda));
    [~,perm] = sort(lambda,'descend');
    lambda = log(lambda(perm).^(1/epsilon));
    psi = psi(:,perm);
    
    % U = S^{-1}Uhat
    psi = Sinv*psi;
    
else
    
    psi = 0;
    lambda = 0;
    
end

%figure(10);plot(-log(diag(b)));hold on;

%%% Normalize qest into a density by dividing by m0
qest = qest/(N*(4*pi*epsilon)^(dim/2));

%%% U^TU = S^{-2} implies that U^TUS^2 = I
%%% componentwise, this means \sum_j phi_j(x)^2 S^2(i,i) = 1 = <phi_j,\phi_j>_peq
%%% but since phi_j is evaluated at x_i with sampling measure q(x),
%%% then S^2(i,i) = peq(x_i)/q(x_i) which means, peq = S^2*q

peqoversample = Ssquare;
peq = qest.*Ssquare;          %%% Invariant measure of the system

%%% normalization factor Z = \frac{1}{N}sum_{j=1}^N peq(x_j)/qest(x_j)
peq = peq./mean(peq./qest);         %%% normalization factor

%%% normalize eigenfunctions such that \sum_i psi(x_i)^2 p(x_i)/q(x_i) = 1
%%% this gives
for i = 1:nvars
    psi(:,i) = psi(:,i)/sqrt(mean(psi(:,i).^2.*(peq./qest)));
end


end




