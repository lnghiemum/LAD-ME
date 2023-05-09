%%% The code implements the likelihood-based estimators for the central
%%% subspace when covariates are contaminated by additive measurement
%%% errors. 
%%% References: Nghiem et al. (2023), Likelihood-based surrogate dimension
%%% reduction

%% Adding paths to current directory
addpath(genpath("./tools"))
addpath(genpath("./manopt"))
addpath(genpath("./SDR_test_rev03092010")); % From 
clear
%% Generating data from a single index model

n = 500; p = 20; 
% Generate matrix for true covariate X

SigmaX = zeros(p, p);
for i=1:p
for j=1:p
    SigmaX(i, j) = 0.5^(abs(i-j));
end
end
X = mvnrnd(zeros(p, 1), SigmaX, n);

% Generate matrix for measurement errors U
SigmaU = diag(unifrnd(0.05, 0.20, p, 1));
U = mvnrnd(zeros(p, 1), SigmaU, n);

% Generate surrogate matrix W
W = X + U;

% Generate y based on a single index model
truebeta =  [ones(3,1); zeros(p-3, 1)];
y = 0.5 * (X * truebeta).^3 + 0.25*randn(n, 1);

%% Computing estimators
% Assuming the true d is known
d = 1;
% Invariance-law estimators based on the adjusted surrogates
hatL = eye(p) - SigmaU * inv(cov(W));
Xhat = W * hatL';

% Slicing y 
nslices = 5;
type = "cont";

% Inverse-moment-based estimators, i.e IL-SIR
[Gammahat_IL_SIR] = computeSIR(Xhat, y, nslices, type, d);
    
% IL-LAD estimators (unpenalized)
Gammahat_IL_LAD = get_lad(Xhat, y, nslices, d , K, options, []);

%% Corrected LAD estimators (unpenalized)
% Set options for optimization on Grassmann manifold
K = speye(p^2) + Tvectranspose(p, p);
options = struct();
options.verbosity = 0;
options.maxiter = 500;
options.tolgradnorm = 1e-10; 

% unpenalized estimator corresponds to the tuning parameter alpha=0
[Gammahat_cLAD] = scLAD(W, y, SigmaU, nslices, d, 0, K, options, [], []);

%% sparse corrected LAD estimators with tuning parameters selected from projection information criterion
alpha_seq = logspace(-3, 0, 40);    
[~, ~, Gammahat_scLAD] = scLAD(W, y, SigmaU, nslices, d, alpha_seq, K, options, Gammahat_cLAD, []);
    
    
    