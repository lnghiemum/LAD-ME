function fit = get_regress(FF,X)
%
% fit = get_regress(FF,X);
%
% This function regresses data X onto basis FF.
%
% =============================================
n = size(X,1);
p = size(X,2);
fit = zeros(n,p);
ff = [ones(rows(FF),1) FF];
for j=1:p,
    Xcent = X(:,j) - mean(X(:,j));
    fit(:,j) = ff*regress(Xcent,ff);
end
