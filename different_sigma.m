clear;
S=30; % Number of end-users per operator
O=16; % Number of operators
N=S*O; % Number of users
Rtot= 0.5*N; % Total resource
%x2=Rtot/N*ones(N,1); % Initial values of allocations for 2-layer
%x3=Rtot/N*ones(S,O); % Initial values of allocations for 3-layer
%x3dc=Rtot/N*ones(S,O); % Initial used values of allocations for 3-layer decentralized
x3adm=Rtot/N*ones(S,O); % Initial for 3-layer ADMM-like
%x3admp=Rtot/N*ones(S,O); % Initial for 3-layer ADMM-like fastest beta activation
%x3admp_inServer = Rtot/N*ones(S,O); 


lambda=1; % Initial lambda
mu=lambda*ones(O,1); % Initial mu
mudc=mu;
lambdadc=mu;
muadm=mu;
lambdaadm=mu;
muadmp=mu;
lambdaadmp=mu;
y=Rtot/O*ones(O,1); % Initial y
ydc=y; % Initial used y
yadm=y;
yadmp=y;

%Graph
k1=4;
k2=4;
Niter = 2000;
w=rand(N,1)+1; % initial weights
h = sqrt(1/2)*(randn(N,1)+1i*randn(N,1)); % initial channel
a = abs(h.^2);

T=(1)*(Niter); % duty cycle
p=lambda;
    
sigma_values = [0.2, 0.4, 0.6, 0.8];
fadm_sigma  = zeros(length(sigma_values), k1*k2*Niter);
%fadm_sigma_err  = zeros(length(sigma_values), k1*k2*Niter);
lat_adm_sigma  = zeros(length(sigma_values), k1*k2*Niter);
algConn = zeros(length(sigma_values),1);

for ind_sigma = 1:length(sigma_values)
    sigma = sigma_values(ind_sigma); 
    prob_of_Connection=sigma;
    W=erdos_renyi(O,prob_of_Connection);
    adjacency = double (W > 0);
    idx = W ~=0;
    c = sum(idx,1);
    max_Neighbor = max(c)-1;
    min_Neighbor = min(c)-1;
    while min_Neighbor == 0
        W=erdos_renyi(O,prob_of_Connection);
        adjacency = double (W > 0);
        idx = W ~=0;
        c = sum(idx,1);
        max_Neighbor = max(c)-1;
        min_Neighbor = min(c)-1;
    end
    G= graph(W,'omitselfloops');
    plot(G)  % manually check if the graph is strongly connected
    edges=table2array(G.Edges);

    %Edge variables
    theta = zeros(size(edges,1),1);
    thetap = zeros(size(edges,1),1);


    % number of iterations and step sizes
    %k1=2;
    %k2=2;
    g2=2/105;
    g1=10;  
    step2=5/100;
    ds_dm=1/5;
    phi=1;
    stepTheta = 1;
    stepLambdaadm = 1/100;

    %Niter=200;

   


    %ADMM random activation probabilities or portion of fastest agents
    %p_edge = 1; % probability of activation
    %p_server = 0.6;
    beta_device = 0.7; % portion of fastest agents
    beta_server = 0.7;

    avg_latency_device = 1;
    avg_latency_server = 100;


    %Save vectors
    x2save=zeros(N,Niter);
    ysave=zeros(O,Niter*k1);
    ydcsave=zeros(O,Niter*k1);
    yadmsave=zeros(O,Niter*k1);
    yadmpsave=zeros(O,Niter*k1);
    ytemp=zeros(O,k1);
    ydctemp=zeros(O,k1);
    musave=zeros(O,Niter*k1*k2);
    mudcsave=zeros(O,Niter*k1*k2);
    muadmsave=zeros(O,Niter*k1*k2);
    muadmpsave=zeros(O,Niter*k1*k2);
    mutemp2=zeros(O,k1*k2);
    mudctemp2=zeros(O,k1*k2);
    mutemp=zeros(O,k2);
    mudctemp=zeros(O,k2);
    psave=zeros(1,Niter);
    lambdasave=zeros(1,Niter);
    lambdadcsave=zeros(O,Niter);
    lambdadmsave=zeros(O,Niter);
    lambdadmpsave=zeros(O,Niter);


    thetasave=zeros(length(theta),Niter);
    thetapsave=zeros(length(theta),Niter);

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

%tol=10
%ptemp=0;

%%

    for i=1:Niter-1
        %% Communication Times
        device_times = exprnd(avg_latency_device,S,O,k1*k2);
        %server_times = exprnd(avg_latency_server,size(edges,1),1);
        server_times = wblrnd(avg_latency_server/2,0.5,size(edges,1),1);
        %server_times = lognrnd(10,1,size(edges,1),1);
        %server_times = rand(size(edges,1),1) * 5 + avg_latency_server;
        
        w3=reshape(w,[S,O]);
        a3=reshape(a,[S,O]);
        x2=w./p - 1./a;
        x2=min(max(0,x2),Rtot);
    
        %% ADMM like method
    
        lat_server = max(server_times);
        lat_adm_server(i) = lat_server;


        for e = 1:size(edges,1)
            theta(e) = theta(e) + stepTheta * (lambdaadm(edges(e,2))-lambdaadm(edges(e,1)));
        end
        thetasave(:,i) = theta;


        for o=1:O
            thetaind1 = find(edges(:,1) == o);
            thetaind2 = find(edges(:,2) == o);
            sumtheta = sum(theta(thetaind1)) - sum(theta(thetaind2));
            lambdaadm(o)=lambdaadm(o)+stepLambdaadm*(yadm(o)-Rtot/O+sumtheta);
            lambdaadm(o)=max(0,lambdaadm(o));

            lambdadmsave(o,i)=lambdaadm(o);

            for j=1:k1
                yadm(o)=yadm(o)+g1*(muadm(o)-lambdaadm(o));
                yadm(o)=max(0,yadm(o));
                ytemp(o,j)=yadm(o);

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

                    muadm(o)=muadm(o)+g2*(sum(x3adm(:,o))-yadm(o));
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
        yadmsave(:,((i-1)*k1+1:(i*k1))) = ytemp;
        muadmsave(:,((i-1)*k1*k2+1:(i*k1*k2))) = mutemp2;
    end

    %f3=sum(f3oj);
    %f3converge=sum(f3cj);
    fadm=sum(f3admoj);
    %fadmp=sum(f3admpoj);
    
    fadm_sigma(ind_sigma,:) = fadm;
    
        %% Run-time and figure
    %lat_admp  = zeros(1,k1*k2*Niter);
    lat_adm  = zeros(1,k1*k2*Niter);
    %lat_allRed = zeros(1,k1*k2*Niter);
    %cur_lat_admp = 0;
    cur_lat_adm = 0;
    %cur_lat_AR = 0;
    %cur_lat_admp_perserver = zeros(O,1);
    cur_lat_adm_perserver = zeros(O,1);
    %cur_lat_AR_perserver = zeros(O,1);

    for ind_lat = 1:k1*k2*Niter
        if mod(ind_lat,k1*k2) == 0
            cur_lat_adm = max(cur_lat_adm_perserver);
            %cur_lat_admp = max(cur_lat_admp_perserver);
            %cur_lat_AR = max(cur_lat_AR_perserver);
            %cur_lat_admp = cur_lat_admp + lat_admp_server(ind_lat/(k1*k2));
            cur_lat_adm = cur_lat_adm + lat_adm_server(ind_lat/(k1*k2));
            %cur_lat_AR = cur_lat_AR + lat_AllReduce_server(ind_lat/(k1*k2));
            %cur_lat_admp_perserver = ones(O,1) * cur_lat_admp;
            cur_lat_adm_perserver = ones(O,1) * cur_lat_adm;
            %cur_lat_AR_perserver = ones(O,1) * cur_lat_AR;
        end
        %cur_lat_admp_perserver = cur_lat_admp_perserver +lat_admp_device(:,ind_lat);
        %lat_admp(ind_lat) = mean(cur_lat_admp_perserver);
        cur_lat_adm_perserver = cur_lat_adm_perserver + lat_adm_device(:,ind_lat);
        lat_adm(ind_lat) = mean(cur_lat_adm_perserver);
        %cur_lat_AR_perserver = cur_lat_AR_perserver + lat_adm_device(:,ind_lat);
        %lat_allRed(ind_lat) = mean(cur_lat_AR_perserver);
    end
    lat_adm_sigma(ind_sigma,:) = lat_adm;
    algConn(ind_sigma) = algebraicConnectivity(adjacency);
    
end

optim = fadm_sigma(length(sigma_values),k1*k2*(Niter-1));
fadm_sigma_err = abs(fadm_sigma - optim); 

% f2longIter=zeros(k1*k2,1);
% f2convergelongIter=zeros(k1*k2,1);
% for i=1:Niter
%     f2longIter(:,i)=f2(i)*ones(k1*k2,1);
%     f2convergelongIter(:,i)=f2converge(i)*ones(k1*k2,1);
% end
% f2longIter=reshape(f2longIter,[1,k1*k2*Niter]);
% f2convergelongIter=reshape(f2convergelongIter,[1,k1*k2*Niter]);
%% Conv Figs
figure;
%plot((1:length(f3converge))/(k1*k2),f3converge,(1:k1*k2*Niter)/(k1*k2),f2convergelongIter,...
%    (1:length(f3))/(k1*k2),f3,(1:length(fadm))/(k1*k2),fadm);
plot(fadm_sigma');

figure;
subplot(1,2,1);
semilogy(fadm_sigma_err');
subplot(1,2,2);
for ind_sigma = 1:length(sigma_values)
    semilogy(lat_adm_sigma(ind_sigma,:), fadm_sigma_err(ind_sigma,:));
    hold on
end

%opt=f2convergelongIter(Niter-1);
%figure;
%semilogy((1:length(f3converge))/(k1*k2),abs(f3converge-opt),(1:k1*k2*Niter)/(k1*k2),abs(f2convergelongIter-opt),...
   % (1:length(f3))/(k1*k2),abs(f3-opt),(1:length(fadm))/(k1*k2),abs(fadm-opt),(1:length(fadmp))/(k1*k2),abs(fadmp-opt));




% figure;
% semilogy(lat_allRed/1000,abs(f3converge-opt)/opt,lat_adm/1000,abs(f3-opt)/opt,lat_adm/1000,abs(fadm-opt)/opt,lat_admp/1000,abs(fadmp-opt)/opt);
% legend('Ring AllReduce','Metropolis','Dec. ML-NUM',strcat('Heur. ','\beta_e= ',num2str(beta_server),',\beta_d= ',num2str(beta_device)));
% 
% figure;
% semilogy((1:length(f3converge))/(k1*k2),abs(f3converge-opt)/opt,(1:length(f3))/(k1*k2),abs(f3-opt)/opt,...
%     (1:length(fadm))/(k1*k2),abs(fadm-opt)/opt,(1:length(fadmp))/(k1*k2),abs(fadmp-opt)/opt);
% legend('Ring AllReduce','Metropolis','Dec. ML-NUM',strcat('Heur. ','\beta_e= ',num2str(beta_server),',\beta_d= ',num2str(beta_device)));
% 
