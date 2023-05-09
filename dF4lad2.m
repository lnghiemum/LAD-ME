function [df] = dF4lad2(W,FParameters, alpha, Tvectranspose)
%	Derivative of F (minus the log-likelihood) for the LAD model.
% Inputs:
%    - W: orthogonal basis matrix for the dimension reduction subspace.
%    - FParameters: structure of parameters computed from the sample. It
%    contains:
%          - FParameters.sigma = array of conditional covariance matrices
%          - FParameters.sigmag = marginal covariance matrix
%          - FParameters.n: sample size for each value of teh response Y.
%   - alpha: tuning parameters
%   - K: a matrix K appears in the gradient computation of penalty term
%
%==========================================================================
    sigma = FParameters.sigma;
    sigmag = FParameters.sigmag;
    p = cols(sigma);
    nj = FParameters.n;
    n = sum(nj);
    f = nj/n;
    a = zeros(rows(W), cols(W), length(f));
    sigma_i = zeros(p);
    for i=1:length(f)
        sigma_i(:,:) = sigma(i,:,:);
        a(:,:,i) = -f(i)*sigma_i*W*inv(W'*sigma_i*W);
    end
    first = sum(a,3);
    second = sigmag * W * inv(W' * sigmag * W); % dim: p * d
    
    % ---Derivative of likelihood function for the LAD model
    llh = -(first+second);
    
    % --- Subdifferential with respect to the penalty 
    p1 = sign(reshape(W * W', p^2, 1)); % dim: p^2  \times 1
    p2 = Tvectranspose * kron(W, eye(p)); % dim: p^2 \times pd
    pensub = alpha .* reshape(p2'*p1, rows(W), cols(W)); % dim: p * d
    df = llh + pensub;
end