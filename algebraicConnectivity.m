function a=algebraicConnectivity(adj)

L=diag(sum(adj))-adj;

[~,D]=eig(L);
s=-sort(-diag(D));
a=s(length(s)-1);