function K = Tvectranspose(m, n) 
% K = Tvectranspose(m,n)  
% The function finds the matrix K of dimension mn * mn such that for any matrix A of dimension m * n, 
% vec(A) = K vec(A^\top)
q = m * n; 

ii = repmat(1:q, 1, q);
jj = repelem(1:q, q);

kk = jj == 1 + m*(ii-1) - (q-1)*floor((ii-1)/n);
K = reshape(kk, q, q);

K = sparse(K);
end