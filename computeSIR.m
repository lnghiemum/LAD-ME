function betaSIR = computeSIR(X, y, h, type, d)
% The function compute SIR estimates, input X = predictor, Y = response,
% and h being the number of slices
nrows = size(X, 1);
n = nrows;
ncols = size(X,2);
if nargin < 5
    d = ncols;
end    

sizes = zeros(h, 1);
means = zeros(ncols, h);

nsqrtx = invsqrtm(cov(X)); %Compute inverse sqrt of marginal covariance matrix
Z = (X - kron(ones(n,1),mean(X)))*nsqrtx ;  

Y = slices(y, h, type);
for j=1:h
    Zj=Z(Y==j,:);
    sizes(j) = size(Zj,1);
    means(:,j) = mean(Zj)';
end

mat2 = zeros(ncols, ncols);
for i=1:h
	mat2 = mat2 + means(:,i)*means(:,i)'*sizes(j)/nrows;
end

S1 = firsteigs(mat2,d);
sir = nsqrtx*S1;
betaSIR = orth(sir);
end