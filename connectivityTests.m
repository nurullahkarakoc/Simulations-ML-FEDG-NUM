clear;
N=50;
Niter=1000;
prob_of_Connection=(0.05: 0.05: 1);
succ = zeros(1,length(prob_of_Connection));
for iter = 1:Niter
    for i=1:length(prob_of_Connection)
        W=erdos_renyi(N,prob_of_Connection(i));
        idx = W ~=0;
        G= graph(W,'omitselfloops');
        %plot(G)
        bins=conncomp(G);
        if sum(bins) == N
            succ(i) = succ(i) +1;
        end
    end
end

%figure;
%plot(prob_of_Connection,succ);

%%
c = [ones(1,N/5) zeros(1,4*N/5)];
A = toeplitz(c,c);
GG = graph(A,'omitselfloops');
plot(GG)

succ_toep = zeros(1,length(prob_of_Connection));
for iter = 1:Niter
    for i=1:length(prob_of_Connection)
        W=erdos_renyi_modified(N,prob_of_Connection(i),A);
        idx = W ~=0;
        G= graph(W,'omitselfloops');
        %plot(G)
        bins=conncomp(G);
        if sum(bins) == N
            succ_toep(i) = succ_toep(i) +1;
        end
    end
end

%%
A = bucky;
GGG = graph(A,'omitselfloops');
figure;
plot(GGG)

succ_circ = zeros(1,length(prob_of_Connection));
for iter = 1:Niter
    for i=1:length(prob_of_Connection)
        W=erdos_renyi_modified(60,prob_of_Connection(i),A);
        idx = W ~=0;
        G= graph(W,'omitselfloops');
        %plot(G)
        bins=conncomp(G);
        if sum(bins) == 60
            succ_circ(i) = succ_circ(i) +1;
        end
    end
end

figure;
plot(prob_of_Connection,succ,prob_of_Connection,succ_toep,prob_of_Connection,succ_circ);
