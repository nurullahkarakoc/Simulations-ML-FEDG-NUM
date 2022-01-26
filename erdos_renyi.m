function [W] = erdos_renyi(N,p)
A= double(rand(N,N)<p);
A=A-diag(diag(A));
B=triu(A);
A=B+B';
A=eye(N)+A;
deg=sum(double(A>0),2);
W=zeros(N,N);
for i=1:N
    for j=1:N
        if i~=j && A(i,j)~=0 
        W(i,j)=1/(max(deg(i),deg(j)));
        end
    end
end
Wtemp=W;
sumW=sum(Wtemp);
for i=1:N
    W(i,i)=1-sumW(i);
end