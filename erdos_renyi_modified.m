function [A] = erdos_renyi_modified(N,p,init)
A= double(rand(N,N)<p).*init;
A=A-diag(diag(A));
B=triu(A);
A=B+B';
A=eye(N)+A;