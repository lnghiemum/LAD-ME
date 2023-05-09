function C = projection(A, B)
C = A * linsolve(A'* B* A,  A'* B);
end