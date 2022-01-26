m=4;
L=10;
d_0=100;
k=5;
ep=200;
alpha=linspace(0,0.15,10000);
rate=max(abs(1-alpha.*m),abs(1-alpha.*L));
bound=rate.^k*d_0+alpha.*ep.*(1-rate.^k)./(1-rate);
[~,index]=min(bound);
semilogy(alpha,bound)
alpha(index)

% syms m L d_0 k ep alpha;
% rate=max(abs(1-alpha*m),abs(1-alpha*L))
% bound=rate^k*d_0+alpha.*ep*(1-rate^k)/(1-rate)
% 
% boundDer=diff(bound,alpha)