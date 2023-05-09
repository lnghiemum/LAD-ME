function [Gammahatest, data_parameters_corrected, Gammahat_corrected_PIC, cost, info] = scLAD(W, y, Sigma_u, nslices, numdir, alpha_seq, K, options, Gammahat0, L)
% the function estimates surrogate SDR from LAD method
% Input:
%    - W: covariate surrogate
%    - y: response vector (assuming continuous for this version)
%    - Sigma_u: measurement error covariance
%    - nslices: number of slices, for this version
%    - numdir: number of dimensions
%    - alpha_seq: sequence of tuning parameters
%    - K: a sparse matrix appearing on the derivation of the \ell-1 norm
%    - Gammahat0: initial estimate        
% Output:
%    - Gammahatest: a p \times \numdir \times length(alpha_seq) array

%% Compute the naive LAD estimator without any penalization
[n,p] = size(W);
d = numdir;
[Gammahat_naive, Deltahat_naive, data_parameters_naive] = get_lad(W, y, nslices, d, K, options, Gammahat0);
%[~, Deltahat_naive, data_parameters_naive, GammahatnaivePIC] = slad(W, y, nslices, 1, alpha_seq, K, options, Gammahat_unpen);

% Corrected lad model
if isempty(L)
estDeltahat = Deltahat_naive - Sigma_u;
L = estDeltahat /(estDeltahat + Sigma_u);
end
Gammahat_corrected_MM = linsolve(L', Gammahat_naive);
Gammahat_corrected_MM = orth(Gammahat_corrected_MM);

% Using MM estimator as the input for the optimization problem
data_parameters_corrected = data_parameters_naive; 
data_parameters_corrected.sigmag = L * data_parameters_naive.sigmag * L';
for j =1:nslices
    data_parameters_corrected.sigma(j, :, :) =  L * squeeze(data_parameters_naive.sigma(j, :, :)) * L'; 
end

Gammahatest = zeros(p, numdir, length(alpha_seq));
for k=1:length(alpha_seq)
  %  disp(k)
manifold = grassmannfactory(p, numdir);
problem.M = manifold;
problem.cost = @(x) n^2*F4lad2(x, data_parameters_corrected, alpha_seq(k));
problem.egrad = @(x) n^2*dF4lad2(x, data_parameters_corrected, alpha_seq(k), K);

% Numerically check gradient consistency (optional).
%checkgradient(problem);

% Solve
warning('off', 'manopt:getHessian:approx');
if isempty(Gammahat0)==1
    Gammahat0 = Gammahat_corrected_MM;
end   
[Gammahat, cost, info, options] = trustregions(problem, Gammahat0, options);

%Gammahat(abs(Gammahat) < 1e-3) = 0;
Gammahatest(:, :, k) = Gammahat;
clear problem
end

% using PIC-type criteria to select the tuning parameters
Gammahat_corrected_PIC = [];
if length(alpha_seq) >1
    loss = zeros(length(alpha_seq), 1);
    df = zeros(length(alpha_seq), 1);
    for kk=1:length(alpha_seq)
        loss(kk) =  norm(computeproj(Gammahat0) - computeproj(Gammahatest(:, :, kk )), "fro");
        df(kk) = sum(abs(diag(Gammahatest(:,:, kk) * Gammahatest(:,:, kk)') ) > 1e-7);
    end
    PIC = loss.^2 + log(p)/p * d.*(df-d);
    Gammahat_corrected_PIC = Gammahatest(:, :, PIC == min(PIC));
end
end