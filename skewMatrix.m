A = [0   1  -1;
    -1   0   1;
    1  -1  0 ];
norm(-ones(1,3)*A^2*ones(3,1))
A = -A^2
s = svd(A)
rank(A)
cond(A)