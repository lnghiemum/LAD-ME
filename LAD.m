function [Gammahat, hatDelta, data_parameters, Deltahat_y, cost, info] = LAD(X, y, nslices, numdir, K, options, Gammaest0)
% Function to compute lad, no penalization
[n, p] = size(X);

parameters.nslices = nslices;
if nslices == length(unique(y))
    type = 'disc';
else 
    type = 'cont';
end
Y = slices(y,parameters.nslices, type);

% get sample statistics, including marginal covariance, and covariance
% within each slices
data_parameters = setdatapars_v2(Y,X,parameters.nslices);

Deltatilde = zeros(p);
for j = 1:nslices
    Deltatilde = Deltatilde + data_parameters.n(j)/n * squeeze(data_parameters.sigma(j, :, :));
end

% Specify the number of dimensions, d = 1, in this case
d = numdir; 

estbeta = zeros(p, numdir);

manifold = grassmannfactory(p, d);
problem.M = manifold;
problem.cost = @(x) F4lad2(x, data_parameters, 0);
problem.egrad = @(x) dF4lad2(x, data_parameters, 0, K);

% Numerically check gradient consistency (optional).
%checkgradient(problem);

% Solve.
warning('off', 'manopt:getHessian:approx');
[Gammahat, cost, info] = trustregions(problem, Gammaest0, options);

invDeltahat = inv(data_parameters.sigmag) + ...
		    Gammahat * linsolve(Gammahat'* Deltatilde * Gammahat, Gammahat') - ...
		    Gammahat * linsolve(Gammahat'* data_parameters.sigmag * Gammahat, Gammahat');
%Gammahat(abs(Gammahat)<1e-3) = 0;        

hatDelta = inv(invDeltahat);

% calculating log likelhood
Pj = projection(Gammahat, hatDelta);
Deltahat_y = zeros(p, p, nslices);
for kk = 1:nslices
    Deltahat_y(:, :, kk) = hatDelta + Pj' * (squeeze(data_parameters.sigma(kk, :, :)) - hatDelta) *Pj;
end

   
end