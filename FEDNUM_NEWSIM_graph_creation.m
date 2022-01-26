clear;

N=50; % it should be even number

Ww=erdos_renyi(N,0.05);
idx = Ww ~=0;
c = sum(idx,1);
max_Neighbor = max(c)-1;
min_Neighbor = min(c)-1;
    while min_Neighbor == 0
        Ww=erdos_renyi(N,0.05);
        idx = Ww ~=0;
        c = sum(idx,1);
        max_Neighbor = max(c)-1;
        min_Neighbor = min(c)-1;
    end

  
G= graph(Ww,'omitselfloops');
h=plot(G,'LineWidth',1.5);  % manually check if the graph is strongly connected
highlight(h,(1:N/2));
edges=G.Edges;
edges_arr=table2array(edges);
L=size(edges,1); 
edgeusage=zeros(L,N);


for i = 1:N/2 % create traffic from i to N-i+1. First half are sources, others destination. 
[P,d,edgepath] = shortestpath(G,i,N-i+1,'Method','unweighted');
edgeusage(edgepath,i)=1; 
end
[M,I] = max(sum(edgeusage,2));
most_used_edge = edges_arr(I(1),:);
highlight(h,most_used_edge(1,1),most_used_edge(1,2),'EdgeColor','r','LineWidth',4);
use_sources = find(edgeusage(I(1),:));
highlight(h,use_sources,'NodeColor','r');

highlight(h,(N/2+1:N),'Marker','^','MarkerSize',8,'NodeColor','k');

str1=compose("s%d",(1:N/2));
str2=compose("d%d",(1:N/2));
labelnode(h,(1:N/2),str1);
labelnode(h,(N:-1:N/2+1),str2);
h.NodeFontSize = 10;


I=abs(full(incidence(G)));

A=zeros(L,L);
for i = 1:L
    for j = 1:N
        if I(j,i) == 1 
            k = find(I(j,:));
            A(i,k) = 1;
        end
    end
end
Gedge=graph(A,'omitselfloops');
figure;
plot(Gedge)
edg=table2array(Gedge.Edges);



deg=sum(double(A>0),2);
W=zeros(L,L);

for i=1:L
    for j=1:L
        if i~=j && A(i,j)~=0 
        W(i,j)=1/(max(deg(i),deg(j)));
        end
    end
end
Wtemp=W;
sumW=sum(Wtemp);
for i=1:L
    W(i,i)=1-sumW(i);
end

WW= W-ones(L,L)/L;    
max(abs(eig(WW)))

