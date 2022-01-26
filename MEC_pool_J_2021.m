clear;
S=30; % Number of end-users per operator
O=10; % Number of operators
N=S*O; % Number of users
P=5; % Number of MEC pools
Rtot= 0.5*N; % Total resource
Cap = [10;10;10;60;60]; %pool caps
%Cap = ones(P,1)*Rtot/P;
pool_to_server= [1 1 1 1 0 0 0 0 0 0;
                 0 0 1 1 1 1 0 0 0 0;
                 0 0 0 0 1 1 1 1 0 0;
                 0 0 0 0 0 0 1 1 1 1;
                 1 1 0 0 0 0 0 0 1 1];
%pool_to_server = double(rand(P,O)>0.2); 
             A = [zeros(P,P) pool_to_server;
                 pool_to_server' zeros(O,O)];
             h=plot(graph(A),'LineWidth',1.5);
             h.YData(1:P) = 2;
h.YData((P+1):end) = 1;
h.XData(1:P) = linspace(0,1,P);
h.XData((P+1):end) = linspace(0,1,O);
%pool_to_server=ones(P,O);
highlight(h,1:5,'Marker','s','MarkerSize',16,'NodeColor','#77AC30');
str1=compose("  P%d",(1:5));
labelnode(h,1:5,str1);
h.NodeFontSize = 12;
labelnode(h,6:15,1:10);



%% 
%x2=Rtot/N*ones(N,1); % Initial values of allocations for 2-layer
x3=Rtot/N*ones(S,O); % Initial values of allocations for 3-layer
x3dc=Rtot/N*ones(S,O); % Initial used values of allocations for 3-layer decentralized
x3adm=Rtot/N*ones(S,O); % Initial for 3-layer ADMM-like
x3admp=Rtot/N*ones(S,O); % Initial for 3-layer ADMM-like fastest beta activation
x3admp_inServer = Rtot/N*ones(S,O);

lambda_value=1;
lambda=lambda_value*ones(P,1); % Initial lambda
mu=lambda_value*ones(O,1); % Initial mu
mudc=mu;
lambdadc=lambda_value*ones(O,P);
muadm=mu;
lambdaadm=lambdadc;
muadmp=mu;
lambdaadmp=lambdadc;
y=Rtot/(O*P)*ones(O,P); % Initial y
ydc=y; % Initial used y
yadm=y;
yadmp=y;

% Metropolis Graph
prob_of_Connection=0.6;
W=erdos_renyi(O,prob_of_Connection);
idx = W ~=0;
c = sum(idx,1);
max_Neighbor = max(c)-1;
min_Neighbor = min(c)-1;
while min_Neighbor == 0
    W=erdos_renyi(O,prob_of_Connection);
    idx = W ~=0;
    c = sum(idx,1);
    max_Neighbor = max(c)-1;
    min_Neighbor = min(c)-1;
end

% Alternate Graph
WW=zeros(O,O);
for i = 1:P
    conn_servers = find(pool_to_server(i,:));
    combs = nchoosek(conn_servers,2);
     for j = 1:length(combs)
        WW(combs(j,1),combs(j,2))=1;
        WW(combs(j,2),combs(j,1))=1;
%         a = mod(combs(j,1)+1,10);
%         b = mod(combs(j,2)+1,10);
%         if a == 0
%             a = 10;
%         end
%         if b == 0
%             b = 10;
%         end
%         WW(a,b)=1;
%         WW(b,a)=1;
    end
end

figure;

G= graph(WW,'omitselfloops');
plot(G)  % manually check if the graph is strongly connected
edges=table2array(G.Edges);

%Edge variables
theta = zeros(size(edges,1),P);
thetap = zeros(size(edges,1),P);

% number of iterations and step sizes
k1=2;
k2=2;
g2=1/105;
g1=5;
g1adm=g1;
g2adm=g2;
step2=1/100;
step2dc=1/100;
ds_dm=1/5;
phi=1;
stepTheta = 2;
stepLambdaadm = 6/100;

Niter=5000;

w=rand(N,1)+1; % initial weights
h = sqrt(1/2)*(randn(N,1)+1i*randn(N,1)); % initial channel
%h=ones(N,1);
a = abs(h.^2);

T=(1)*(Niter); % duty cycle
p=lambda_value;


%ADMM random activation probabilities or portion of fastest agents
%p_edge = 1; % probability of activation
%p_server = 0.6;
beta_device = 0.7; % portion of fastest agents
beta_server = 0.7;

avg_latency_device = 10;
avg_latency_server = 100;


%Save vectors
%x2save=zeros(N,Niter);
ysave=zeros(O,P,Niter*k1);
ydcsave=zeros(O,P,Niter*k1);
yadmsave=zeros(O,P,Niter*k1);
yadmpsave=zeros(O,P,Niter*k1);
ytemp=zeros(O,P,k1);
ydctemp=zeros(O,P,k1);
musave=zeros(O,Niter*k1*k2);
mudcsave=zeros(O,Niter*k1*k2);
muadmsave=zeros(O,Niter*k1*k2);
muadmpsave=zeros(O,Niter*k1*k2);
mutemp2=zeros(O,k1*k2);
mudctemp2=zeros(O,k1*k2);
mutemp=zeros(O,k2);
mudctemp=zeros(O,k2);
%psave=zeros(1,Niter);
lambdasave=zeros(P,Niter);
lambdadcsave=zeros(O,P,Niter);
lambdadmsave=zeros(O,P,Niter);
lambdadmpsave=zeros(O,P,Niter);

[sz1,sz2]=size(theta);
thetasave=zeros(sz1,sz2,Niter);
thetapsave=zeros(sz1,sz2,Niter);

f2=zeros(1,Niter);
f2converge=zeros(1,Niter);
fadm=zeros(1,Niter);
fadmp=zeros(1,Niter);
Res=zeros(1,Niter);
x2u=Rtot/N*ones(N,1);
x3u=Rtot/N*ones(S,O);
f3o=zeros(O,1);
f3oj=zeros(O,k1*k2*Niter);
f3admo=zeros(O,1);
f3admoj=zeros(O,k1*k2*Niter);
f3admpo=zeros(O,1);
f3admpoj=zeros(O,k1*k2*Niter);
f3c=zeros(O,1);
f3cj=zeros(O,k1*k2*Niter);

lat_admp_device = zeros(O,k1*k2*Niter);
lat_admp_server = zeros(1,Niter);

lat_adm_device = zeros(O,k1*k2*Niter);
lat_adm_server = zeros(1,Niter);

lat_AllReduce_server = zeros(1,Niter);
lat_AllReduce_device = zeros(O,k1*k2*Niter);

cap_vec = ones(O,1)*Cap'./O;
% for o =1:O
%     for ind_pts = 1:P
%         if pool_to_server(ind_pts,o) == 0
%             cap_vec(o,ind_pts)=0;
%         end
%     end
% end
%tol=10
%ptemp=0;

%%

for i=1:Niter-1
    %% Communication Times
    device_times = exprnd(avg_latency_device,S,O,k1*k2);
    server_times = exprnd(avg_latency_server,size(edges,1),1);
    %server_times = wblrnd(avg_latency_server/2,0.5,size(edges,1),1);
    %server_times = lognrnd(10,1,size(edges,1),1);
    %server_times = rand(size(edges,1),1) * 5 + avg_latency_server;
    
    
    
    
    %% 2-Layer
    %flag=0;
    %     if(mod(i,T)==0) % renew parameters each time slot
    %         w=rand(N,1)+1;
    %         h = sqrt(1/2)*(randn(N,1)+1i*randn(N,1));
    %         %a = abs(h.^2);
    %         a=ones(N,1);
    %         flag=1;
    %     end
    %     if(abs(p-ptemp)<tol)
    %         flag=1;
    %
    %     end
    w3=reshape(w,[S,O]);
    a3=reshape(a,[S,O]);
         x2=w./p - 1./a;
         x2=min(max(0,x2),Rtot);
    
    %     if sum(y)-Rtot > 0 % used rates
    %         yu=y-((sum(y)-Rtot)/O).*ones(O,1);
    %         yu=min(max(0,yu),Rtot);
    %     else
    %         yu=y;
    %     end
    %  x2u=x2;
    %     if sum(x2)-Rtot > 0
    %         x2u=x2-((sum(x2)-Rtot)/N).*ones(N,1);
    %         x2u=min(max(0,x2u),Rtot);
    %         while sum(x2u)-Rtot > 0
    %         x2u=x2u-((sum(x2u)-Rtot)/N).*ones(N,1);
    %         x2u=min(max(0,x2u),Rtot);
    %         end
    %     else
    %         x2u=x2;
    %     end
    
        x2save(:,i)=x2;
        ptemp=p;
        p=p+step2*(sum(x2)-Rtot);
        p=max(0,p);
        psave(i)=p;
    
        Res(i)=sum(x2);
        % Resu(i)=sum(x2u);
        f2(i)=sum(w.*log(1+a.*x2u));
        f2converge(i)=sum(w.*log(1+a.*x2));
    
    
    
    
    
    %% 3-Layer
    
    lat_server = 2*sum(server_times(1:(O-1),1));
    lat_AllReduce_server(i) = lat_server;
    
    lambda=lambda+step2*(sum(y)'-Cap);
    lambda=max(0,lambda);
    lambdasave(:,i)=lambda;
    
    for o=1:O
        for j=1:k1
            y(o,:)=y(o,:)+g1*(mu(o)*ones(1,P)-lambda');
            y(o,:)=max(0,y(o,:));
            for ind_pts = 1:P
                if pool_to_server(ind_pts,o) == 0
                    y(o,ind_pts)=0;
                end
            end
            ytemp(o,:,j)=y(o,:);
            
            %             if (j==2 && (flag==1 || i==1))
            %                 x3u(:,o)=x3(:,o);
            %             end
            for h=1:k2
                
                x3(:,o)=w3(:,o)./mu(o)-1./a3(:,o);
                x3(:,o)=min(max(0,x3(:,o)),Rtot);
                
                %                 if sum(x3(:,o))-yu(o) > 0
                %                     x3u(:,o)=x3(:,o)-((sum(x3(:,o))-yu(o))/S).*ones(S,1);
                %                     x3u(:,o)=min(max(0,x3u(:,o)),Rtot);
                %                 else
                %                     x3u(:,o)=x3(:,o);
                %                 end
                
                mu(o)=mu(o)+g2*(sum(x3(:,o))-sum(y(o,:)));
                mu(o)=max(0,mu(o));
                % f3o(o)=sum(w3(:,o).*log(1+a3(:,o).*x3u(:,o)));
                f3c(o)=sum(w3(:,o).*log(1+a3(:,o).*x3(:,o)));
                % f3oj(o,(i-1)*k1*k2+(j-1)*k2+h)=f3o(o);
                f3cj(o,(i-1)*k1*k2+(j-1)*k2+h)=f3c(o);
                mutemp(o,h)=mu(o);
                
            end
            mutemp2(o,((j-1)*k2+1:(j*k2))) = mutemp(o,:);
        end
    end
    ysave(:,:,((i-1)*k1+1:(i*k1))) = ytemp;
    musave(:,((i-1)*k1*k2+1:(i*k1*k2))) = mutemp2;
    
    %% 3-Layer DC
    
    
    % Write consensus here
    lambdadc=lambdadc+(step2dc).*(O.*ydc-ones(O,1)*Cap');
    lambdadc=max(0,lambdadc);
    
    
    lambdadc=(W^phi)*lambdadc;
    lambdadcsave(:,:,i)=lambdadc;
    for o=1:O
        for j=1:k1
            ydc(o,:)=ydc(o,:)+g1*(mudc(o)*ones(1,P)-lambdadc(o,:));
            ydc(o,:)=max(0,ydc(o,:));
            for ind_pts = 1:P
                if pool_to_server(ind_pts,o) == 0
                    ydc(o,ind_pts)=0;
                end
            end
            ydctemp(o,:,j)=ydc(o,:);
            
            %             if (j==2 && (flag==1 || i==1))
            %                 x3u(:,o)=x3(:,o);
            %             end
            for h=1:k2
                
                x3dc(:,o)=w3(:,o)./mudc(o)-1./a3(:,o);
                x3dc(:,o)=min(max(0,x3dc(:,o)),Rtot);
                
                %                 if sum(x3(:,o))-yu(o) > 0
                %                     x3u(:,o)=x3(:,o)-((sum(x3(:,o))-yu(o))/S).*ones(S,1);
                %                     x3u(:,o)=min(max(0,x3u(:,o)),Rtot);
                %                 else
                %                     x3u(:,o)=x3(:,o);
                %                 end
                
                mudc(o)=mudc(o)+g2*(sum(x3dc(:,o))-sum(ydc(o,:)));
                mudc(o)=max(0,mudc(o));
                f3o(o)=sum(w3(:,o).*log(1+a3(:,o).*x3dc(:,o)));
                %  f3c(o)=sum(w3(:,o).*log(1+a3(:,o).*x3(:,o)));
                f3oj(o,(i-1)*k1*k2+(j-1)*k2+h)=f3o(o);
                % f3cj(o,(i-1)*k1*k2+(j-1)*k2+h)=f3c(o);
                mudctemp(o,h)=mudc(o);
                
            end
            mudctemp2(o,((j-1)*k2+1:(j*k2))) = mudctemp(o,:);
        end
    end
    ydcsave(:,:,((i-1)*k1+1:(i*k1))) = ydctemp;
    mudcsave(:,((i-1)*k1*k2+1:(i*k1*k2))) = mudctemp2;
    
    %% ADMM like method
    
    lat_server = max(server_times);
    lat_adm_server(i) = lat_server;
    
    
    for e = 1:size(edges,1)
        theta(e,:) = theta(e,:) + stepTheta .*(lambdaadm(edges(e,2),:)-lambdaadm(edges(e,1),:));
    end
    thetasave(:,:,i) = theta;
    
    
    for o=1:O
        thetaind1 = find(edges(:,1) == o);
        thetaind2 = find(edges(:,2) == o);
        sumtheta = sum(theta(thetaind1,:),1) - sum(theta(thetaind2,:),1);
        lambdaadm(o,:)=lambdaadm(o,:)+stepLambdaadm*(yadm(o,:)-Cap'./O +sumtheta);
        lambdaadm(o,:)=max(0,lambdaadm(o,:));
        
        lambdadmsave(o,:,i)=lambdaadm(o,:);
        
        for j=1:k1
            yadm(o,:)=yadm(o,:)+g1adm*(muadm(o)*ones(1,P)-lambdaadm(o,:));
            yadm(o,:)=max(0,yadm(o,:));
            for ind_pts = 1:P
                if pool_to_server(ind_pts,o) == 0
                    yadm(o,ind_pts)=0;
                end
            end
            ytemp(o,:,j)=yadm(o,:);
            
            %             if (j==2 && (flag==1 || i==1))
            %                 x3u(:,o)=x3(:,o);
            %             end
            for h=1:k2
                
                lat_device = max(device_times(:,o,(j-1)*k1+h));
                
                x3adm(:,o)=w3(:,o)./muadm(o)-1./a3(:,o);
                x3adm(:,o)=min(max(0,x3adm(:,o)),Rtot);
                
                %                 if sum(x3(:,o))-yu(o) > 0
                %                     x3u(:,o)=x3(:,o)-((sum(x3(:,o))-yu(o))/S).*ones(S,1);
                %                     x3u(:,o)=min(max(0,x3u(:,o)),Rtot);
                %                 else
                %                     x3u(:,o)=x3(:,o);
                %                 end
                
                muadm(o)=muadm(o)+g2adm*(sum(x3adm(:,o))-sum(yadm(o,:)));
                muadm(o)=max(0,muadm(o));
                % f3o(o)=sum(w3(:,o).*log(1+a3(:,o).*x3u(:,o)));
                f3admo(o)=sum(w3(:,o).*log(1+a3(:,o).*x3adm(:,o)));
                % f3oj(o,(i-1)*k1*k2+(j-1)*k2+h)=f3o(o);
                f3admoj(o,(i-1)*k1*k2+(j-1)*k2+h)=f3admo(o);
                mutemp(o,h)=muadm(o);
                
                lat_adm_device(o,(i-1)*k1*k2+(j-1)*k2+h)=lat_device;
                
            end
            mutemp2(o,((j-1)*k2+1:(j*k2))) = mutemp(o,:);
        end
    end
    yadmsave(:,:,((i-1)*k1+1:(i*k1))) = ytemp;
    muadmsave(:,((i-1)*k1*k2+1:(i*k1*k2))) = mutemp2;
    
    %% ADMM like with random activation
    % Active servers
    number_act_links = ceil(beta_server*size(edges,1));
    sorted_server_times = sort(server_times);
    lat_server = sorted_server_times(number_act_links);
    is_link_active = logical (server_times <= lat_server);
    
    lat_admp_server(i) = lat_server;
    
    for e = 1:size(edges,1)
        %        if rand < p_server
        if is_link_active(e)
            thetap(e,:) = thetap(e,:) + stepTheta .* (lambdaadmp(edges(e,2),:)-lambdaadmp(edges(e,1),:));
        end
    end
    thetapsave(:,:,i) = thetap;
    
    
    for o=1:O
        thetapind1 = find(edges(:,1) == o);
        thetapind2 = find(edges(:,2) == o);
        sumthetap = sum(thetap(thetapind1,:),1) - sum(thetap(thetapind2,:),1);
        lambdaadmp(o,:)=lambdaadmp(o,:)+stepLambdaadm*(yadmp(o,:)-Cap'/O+sumthetap);
        lambdaadmp(o,:)=max(0,lambdaadmp(o,:));
        
        lambdadmpsave(o,:,i)=lambdaadmp(o,:);
        
        for j=1:k1
            yadmp(o,:)=yadmp(o,:)+g1*(muadmp(o)*ones(1,P)-lambdaadmp(o,:));
            yadmp(o,:)=max(0,yadmp(o,:));
            for ind_pts = 1:P
                if pool_to_server(ind_pts,o) == 0
                    yadmp(o,ind_pts)=0;
                end
            end
            ytemp(o,:,j)=yadmp(o,:);
            
            for h=1:k2
                
                % Active devices
                number_act_devices = ceil(beta_device*S);
                sorted_device_times = sort(device_times(:,o,(j-1)*k1+h));
                lat_device = sorted_device_times(number_act_devices);
                is_device_active = logical (device_times(:,o,(j-1)*k1+h) <= lat_device);
                
                x3admp(:,o)=w3(:,o)./muadmp(o)-1./a3(:,o);
                x3admp(:,o)=min(max(0,x3admp(:,o)),Rtot);
                %if rand < p_edge
                for s=1:S
                    if is_device_active(s)
                        x3admp_inServer(s,o) = x3admp(s,o);
                    end
                end
                
                
                muadmp(o)=muadmp(o)+g2*(sum(x3admp_inServer(:,o))-sum(yadmp(o,:)));
                muadmp(o)=max(0,muadmp(o));
                f3admpo(o)=sum(w3(:,o).*log(1+a3(:,o).*x3admp(:,o)));
                f3admpoj(o,(i-1)*k1*k2+(j-1)*k2+h)=f3admpo(o);
                mutemp(o,h)=muadmp(o);
                
                lat_admp_device(o,(i-1)*k1*k2+(j-1)*k2+h)=lat_device;
                
                
            end
            mutemp2(o,((j-1)*k2+1:(j*k2))) = mutemp(o,:);
        end
    end
    yadmpsave(:,:,((i-1)*k1+1:(i*k1))) = ytemp;
    muadmpsave(:,((i-1)*k1*k2+1:(i*k1*k2))) = mutemp2;
end

f3=sum(f3oj);
f3converge=sum(f3cj);
fadm=sum(f3admoj);
fadmp=sum(f3admpoj);

f2longIter=zeros(k1*k2,1);
f2convergelongIter=zeros(k1*k2,1);
for i=1:Niter
    f2longIter(:,i)=f2(i)*ones(k1*k2,1);
    f2convergelongIter(:,i)=f2converge(i)*ones(k1*k2,1);
end
f2longIter=reshape(f2longIter,[1,k1*k2*Niter]);
f2convergelongIter=reshape(f2convergelongIter,[1,k1*k2*Niter]);
%% Conv Figs
figure;
plot((1:length(f3converge))/(k1*k2),f3converge,...
    (1:length(f3))/(k1*k2),f3,(1:length(fadm))/(k1*k2),fadm);
legend('AllRed','Metropolis','FEDg-NUM');

opt=f3converge(Niter*k1*k2-k1*k2);
%opt=fadm(Niter-k1*k2);
%opt=f2convergelongIter(Niter-1);

figure;
semilogy((1:length(f3converge))/(k1*k2),abs(f3converge-opt)/opt,...
    (1:length(f3))/(k1*k2),abs(f3-opt)/opt,(1:length(fadm))/(k1*k2),abs(fadm-opt)/opt)
%(1:length(fadmp))/(k1*k2),abs(fadmp-opt)/opt);
legend('AllRed','Metropolis','FEDg-NUM');

%% Run-time and figure
lat_admp  = zeros(1,k1*k2*Niter);
lat_adm  = zeros(1,k1*k2*Niter);
lat_allRed = zeros(1,k1*k2*Niter);
cur_lat_admp = 0;
cur_lat_adm = 0;
cur_lat_AR = 0;
cur_lat_admp_perserver = zeros(O,1);
cur_lat_adm_perserver = zeros(O,1);
cur_lat_AR_perserver = zeros(O,1);

for ind_lat = 1:k1*k2*Niter
    if mod(ind_lat,k1*k2) == 0
        cur_lat_adm = max(cur_lat_adm_perserver);
        cur_lat_admp = max(cur_lat_admp_perserver);
        cur_lat_AR = max(cur_lat_AR_perserver);
        cur_lat_admp = cur_lat_admp + lat_admp_server(ind_lat/(k1*k2));
        cur_lat_adm = cur_lat_adm + lat_adm_server(ind_lat/(k1*k2));
        cur_lat_AR = cur_lat_AR + lat_AllReduce_server(ind_lat/(k1*k2));
        cur_lat_admp_perserver = ones(O,1) * cur_lat_admp;
        cur_lat_adm_perserver = ones(O,1) * cur_lat_adm;
        cur_lat_AR_perserver = ones(O,1) * cur_lat_AR;
    end
    cur_lat_admp_perserver = cur_lat_admp_perserver +lat_admp_device(:,ind_lat);
    lat_admp(ind_lat) = mean(cur_lat_admp_perserver);
    cur_lat_adm_perserver = cur_lat_adm_perserver + lat_adm_device(:,ind_lat);
    lat_adm(ind_lat) = mean(cur_lat_adm_perserver);
    cur_lat_AR_perserver = cur_lat_AR_perserver + lat_adm_device(:,ind_lat);
    lat_allRed(ind_lat) = mean(cur_lat_AR_perserver);
end

fadmp(Niter*k1*k2-k1*k2+1:Niter*k1*k2) = fadmp(Niter*k1*k2-k1*k2);
fadm(Niter*k1*k2-k1*k2+1:Niter*k1*k2) = fadm(Niter*k1*k2-k1*k2);
f3converge(Niter*k1*k2-k1*k2+1:Niter*k1*k2) = f3converge(Niter*k1*k2-k1*k2);
f3(Niter*k1*k2-k1*k2+1:Niter*k1*k2) = f3(Niter*k1*k2-k1*k2);


figure;
semilogy(lat_allRed/1000,abs(f3converge-opt)/opt,lat_adm/1000,abs(f3-opt)/opt,lat_adm/1000,abs(fadm-opt)/opt);
%legend(strcat('\beta_e = ',num2str(beta_server),'\beta_d = ',num2str(beta_device)),'Dec. ML-NUM','AllRed');
legend('AllRed','Metropolis','FEDg-NUM');

%%
idx = W ~=0;
c = sum(idx,1);
max_Neighbor = max(c)-1;
scale1= N/(O+k1*k2*S*ds_dm);
scale2= N/(max_Neighbor+k1*k2*S*ds_dm);

figure;
subplot(2,1,1);
h=plot(G);
h.NodeFontSize = 12;
subplot(2,1,2);
semilogy(lat_allRed/1000,abs(f3converge-opt)/opt,lat_adm/1000,abs(f3-opt)/opt,lat_adm/1000,abs(fadm-opt)/opt);
%legend(strcat('\beta_e = ',num2str(beta_server),'\beta_d = ',num2str(beta_device)),'Dec. ML-NUM','AllRed');
legend('AllRed','Metropolis','FEDg-NUM');

%figure;
%plot((1:length(f3converge))/(k1*k2*scale1),f3converge,(1:k1*k2*Niter)/(k1*k2),f2convergelongIter,(1:length(f3))/(k1*k2*scale2),f3);

%figure;
%plot((1:Niter),lambdadcsave)

% figure;
% semilogy((1:99),abs(ysave(:,1:99)-ysave(:,99))./ysave(:,99));
% xlabel('Iterations $j$');
% ylabel('$|y_o - y_o^*|/|y_o^*|$');
%
% figure;
% semilogy([1:99],abs(musave(:,1:99)-musave(:,99))./musave(:,99));
% xlabel('Iterations $j$');
% ylabel('$|y_o - y_o^*|/|y_o^*|$');
%
% figure;
% semilogy([1:99],abs(f2converge(:,1:99)-f2converge(:,99))./f2converge(:,99));
% xlabel('Iterations $j$');
% ylabel('$|y_o - y_o^*|/|y_o^*|$');