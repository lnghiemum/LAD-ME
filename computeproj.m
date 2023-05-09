function proj = computeproj(A)
proj = A * linsolve(A'*A, A');
end