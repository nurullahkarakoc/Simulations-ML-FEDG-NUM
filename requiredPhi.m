clear;
S=20; % Number of end-users per operator
O=10; % Number of operators
N=S*O; % Number of users
Rtot= 0.5*N; % Total resource
x2=Rtot/N*ones(N,1); % Initial values of allocations for 2-layer
x3=Rtot/N*ones(S,O); % Initial values of allocations for 3-layer
x3dc=Rtot/N*ones(S,O); % Initial used values of allocations for 3-layer

lambda=1; % Initial lambda
mu=lambda*ones(O,1); % Initial mu
mudc=mu;
lambdadc=mu;
y=Rtot/O*ones(O,1); % Initial y
ydc=y; % Initial used y

W=erdos_renyi(O,0.8);
idx = W ~=0;
c = sum(idx,1);
max_Neighbor = max(c)-1;
min_Neighbor = min(c)-1;
    while min_Neighbor == 0
        W=erdos_renyi(O,0.25);
        idx = W ~=0;
        c = sum(idx,1);
        max_Neighbor = max(c)-1;
        min_Neighbor = min(c)-1;
    end


WW= W-ones(O,O)/O;    
max(abs(eig(WW)))
%N=10;
M=1;
w=0.001;
alpha=0.01;
C=5;
psi=0.5;

phi=(log(w)-log(4*M*O*(alpha*O*C+w)))/log(psi);