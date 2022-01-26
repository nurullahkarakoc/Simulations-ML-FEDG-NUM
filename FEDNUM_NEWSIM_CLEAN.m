N=50;
Ptot= 0.1*L; % Total resource % roughly 100mW transmit power per link
%x2=Rtot/N*ones(N,1); % Initial values of allocations for 2-layer
N=N/2; % we will only iterate the sources
x3=Ptot/N*ones(N,1); % Initial values of allocations for 3-layer
x3dc=Ptot/N*ones(N,1); % Initial used values of allocations for 3-layer
x3adm=x3dc;

lambda=0.5; % Initial lambda
mu=lambda*ones(L,1); % Initial mu
mudc=mu;
muadm=mu;
lambdadc=mu;
lambdaadm=mu;
y=Ptot/L*ones(L,1); % Initial y
ydc=y; % Initial used y
yadm=y;

%Edge variables
theta = zeros(size(edg,1),1);


k1=2;
k2=2;
g2=4/105;
g1=0.5;
step2=2/1000;
ds_dm=1/5;
phi=1;
stepTheta = 2;
stepLambdaadm = 4/100;

Niter=1000;
k=Niter;

w=rand(N,1)+1; % initial weights
h = sqrt(1/2)*(randn(N,1)+1i*randn(N,1)); % initial channel
%a = abs(h.^2);
a=ones(N,1);
T=(1)*(Niter); % duty cycle
%p=lambda;
%w=ones(N,1);

%x2save=zeros(N,Niter);
ysave=zeros(L,Niter*k1);
ydcsave=zeros(L,Niter*k1);
yadmsave=zeros(L,Niter*k1);

ytemp=zeros(L,k1);
ydctemp=zeros(L,k1);
yadmtemp=zeros(L,k1);
musave=zeros(L,Niter*k1*k2);
mudcsave=zeros(L,Niter*k1*k2);
muadmsave=zeros(L,Niter*k1*k2);
mutemp2=zeros(L,k1*k2);
mudctemp2=zeros(L,k1*k2);
muadmtemp2=zeros(L,k1*k2);
mutemp=zeros(L,k2);
mudctemp=zeros(L,k2);
muadmtemp=zeros(L,k2);
%psave=zeros(1,Niter);
lambdasave=zeros(1,Niter);
lambdadcsave=zeros(L,Niter);
lambdaadmsave=zeros(L,Niter);

thetasave = zeros(length(theta),Niter);


f2=zeros(1,Niter);
f2converge=zeros(1,Niter);
Res=zeros(1,Niter);
%x2u=Ptot/N*ones(N,1);
%x3u=Ptot/N*ones(N,1);
f3o=zeros(L,1);
f3oj=zeros(L,k1*k2*k);
f3c=zeros(L,1);
f3cj=zeros(L,k1*k2*k);

avg_latency_device = 10;
avg_latency_server = 100;

util=[];
utildc=[];
utiladm=[];

lat_admp_device = zeros(L,k1*k2*Niter);
lat_admp_server = zeros(1,Niter);

lat_adm_device = zeros(L,k1*k2*Niter);
lat_adm_server = zeros(1,Niter);

lat_AllReduce_server = zeros(1,Niter);
lat_AllReduce_device = zeros(L,k1*k2*Niter);

%tol=10;
%ptemp=0;
%%
for i = 1: Niter
    %% Communication Times
    device_times = exprnd(avg_latency_device,N,1,k1*k2);
    server_times = exprnd(avg_latency_server,size(edg,1),1);
    %server_times = wblrnd(avg_latency_server/2,0.5,size(edges,1),1);
    %server_times = lognrnd(10,1,size(edges,1),1);
    %server_times = rand(size(edges,1),1) * 5 + avg_latency_server;
    
    
    %% Central Lambda
    lat_server = 2*sum(server_times(1:(L-1),1));
    lat_AllReduce_server(i) = lat_server;
    
    lambda=lambda+step2*(sum(y)-Ptot);
    lambda=max(0,lambda);
    lambdasave(i)=lambda;
    
    for j=1:k1
        for o=1:L
            y(o)=y(o)+g1*(mu(o)/(1+y(o))-lambda);
            y(o)=max(0,y(o));
            ytemp(o,j)=y(o);
        end
        
        
        for h=1:k2
            
            sumMu=edgeusage(:,1:N)'*mu;
            
            for s=1:N
                x3(s)=w(s)./(sumMu(s))-1./a(s);
                x3(s)=min(max(0,x3(s)),Ptot);
            end
            
            sumX=edgeusage(:,1:N)*x3;
            
            for o = 1:L
                mu(o)=mu(o)+g2*(sumX(o)-log(1+y(o)));
                mu(o)=max(0,mu(o));
                % f3o(o)=sum(w3(:,o).*log(1+a3(:,o).*x3u(:,o)));
                %f3c(o)=sum(w(:,o).*log(1+a3(:,o).*x3(:,o)));
                % f3oj(o,(i-1)*k1*k2+(j-1)*k2+h)=f3o(o);
                % f3cj(o,(i-1)*k1*k2+(j-1)*k2+h)=f3c(o);
                mutemp(o,h)=mu(o);
            end
            
            util = [util sum(w.*log(1+a.*x3))];
            
        end
        for o = 1:L
            mutemp2(o,((j-1)*k2+1:(j*k2))) = mutemp(o,:);
        end
    end
    
    ysave(:,((i-1)*k1+1:(i*k1))) = ytemp;
    musave(:,((i-1)*k1*k2+1:(i*k1*k2))) = mutemp2;
    
    
    
    %% Consensus - Metropolis
    
    lambdadc=lambdadc+step2*(ydc*L-Ptot);
    lambdadc=max(0,lambdadc);
    
    lambdadc=(W^phi)*lambdadc;
    lambdadcsave(:,i)=lambdadc;
    for j=1:k1
        for o=1:L
            ydc(o)=ydc(o)+g1*(mudc(o)/(1+ydc(o))-lambdadc(o));
            ydc(o)=max(0,ydc(o));
            ydctemp(o,j)=ydc(o);
        end
        
        
        for h=1:k2
            
            sumMudc=edgeusage(:,1:N)'*mudc;
            
            for s=1:N
                x3dc(s)=w(s)./(sumMudc(s))-1./a(s);
                x3dc(s)=min(max(0,x3dc(s)),Ptot);
            end
            
            sumXdc=edgeusage(:,1:N)*x3dc;
            
            for o = 1:L
                mudc(o)=mudc(o)+g2*(sumXdc(o)-log(1+ydc(o)));
                mudc(o)=max(0,mudc(o));
                % f3o(o)=sum(w3(:,o).*log(1+a3(:,o).*x3u(:,o)));
                %f3c(o)=sum(w(:,o).*log(1+a3(:,o).*x3(:,o)));
                % f3oj(o,(i-1)*k1*k2+(j-1)*k2+h)=f3o(o);
                % f3cj(o,(i-1)*k1*k2+(j-1)*k2+h)=f3c(o);
                mudctemp(o,h)=mudc(o);
            end
            
            utildc = [utildc sum(w.*log(1+a.*x3dc))];
            
        end
        for o = 1:L
            mudctemp2(o,((j-1)*k2+1:(j*k2))) = mudctemp(o,:);
        end
    end
    
    ydcsave(:,((i-1)*k1+1:(i*k1))) = ydctemp;
    mudcsave(:,((i-1)*k1*k2+1:(i*k1*k2))) = mudctemp2;
    
    
    
    %% Fedg-NUM
    
    lat_server = max(server_times);
    lat_adm_server(i) = lat_server;
    
    for e = 1:size(edg,1)
        theta(e) = theta(e) + stepTheta * (lambdaadm(edg(e,2))-lambdaadm(edg(e,1)));
    end
    thetasave(:,i) = theta;
    
    for o=1:L
        thetaind1 = find(edg(:,1) == o);
        thetaind2 = find(edg(:,2) == o);
        sumtheta = sum(theta(thetaind1)) - sum(theta(thetaind2));
        lambdaadm(o)=lambdaadm(o)+stepLambdaadm*(yadm(o)-Ptot/L+sumtheta);
        lambdaadm(o)=max(0,lambdaadm(o));
        
        lambdaadmsave(o,i)=lambdaadm(o);
    end
    
    for j=1:k1
        for o=1:L
            yadm(o)=yadm(o)+g1*(muadm(o)/(1+yadm(o))-lambdaadm(o));
            yadm(o)=max(0,yadm(o));
            yadmtemp(o,j)=yadm(o);
        end
        
        for h=1:k2
            
            for o=1:L
                activenodes = find(edgeusage(o,:));
                lat_device = max(device_times(activenodes,1,(j-1)*k1+h));
                if ~isempty(lat_device)
                    lat_adm_device(o,(i-1)*k1*k2+(j-1)*k2+h)=lat_device;
                end
            end
            
            
            sumMuadm=edgeusage(:,1:N)'*muadm;
            
            for s=1:N
                x3adm(s)=w(s)./(sumMuadm(s))-1./a(s);
                x3adm(s)=min(max(0,x3adm(s)),Ptot);
            end
            
            sumXadm=edgeusage(:,1:N)*x3adm;
            
            for o = 1:L
                muadm(o)=muadm(o)+g2*(sumXadm(o)-log(1+yadm(o)));
                muadm(o)=max(0,muadm(o));
                % f3o(o)=sum(w3(:,o).*log(1+a3(:,o).*x3u(:,o)));
                %f3c(o)=sum(w(:,o).*log(1+a3(:,o).*x3(:,o)));
                % f3oj(o,(i-1)*k1*k2+(j-1)*k2+h)=f3o(o);
                % f3cj(o,(i-1)*k1*k2+(j-1)*k2+h)=f3c(o);
                muadmtemp(o,h)=muadm(o);
            end
            
            utiladm = [utiladm sum(w.*log(1+a.*x3adm))];
            
        end
        for o = 1:L
            muadmtemp2(o,((j-1)*k2+1:(j*k2))) = muadmtemp(o,:);
        end
    end
    
    yadmsave(:,((i-1)*k1+1:(i*k1))) = yadmtemp;
    muadmsave(:,((i-1)*k1*k2+1:(i*k1*k2))) = muadmtemp2;
    
    
end
plot(1:length(util),util,1:length(utildc),utildc,1:length(utiladm),utiladm)
opt = util(end);
figure;
semilogy(1:length(util),abs(util-opt)/opt,1:length(utildc),abs(utildc-opt)/opt,1:length(utiladm),...
    abs(utiladm-opt)/opt)

%% Run-time and figure
lat_adm  = zeros(1,k1*k2*Niter);
lat_allRed = zeros(1,k1*k2*Niter);
cur_lat_adm = 0;
cur_lat_AR = 0;
cur_lat_adm_perserver = zeros(L,1);
cur_lat_AR_perserver = zeros(L,1);

for ind_lat = 1:k1*k2*Niter
    if mod(ind_lat,k1*k2) == 0
        cur_lat_adm = max(cur_lat_adm_perserver);
        cur_lat_AR = max(cur_lat_AR_perserver);
        cur_lat_adm = cur_lat_adm + lat_adm_server(ind_lat/(k1*k2));
        cur_lat_AR = cur_lat_AR + lat_AllReduce_server(ind_lat/(k1*k2));
        cur_lat_adm_perserver = ones(L,1) * cur_lat_adm;
        cur_lat_AR_perserver = ones(L,1) * cur_lat_AR;
    end
    
    cur_lat_adm_perserver = cur_lat_adm_perserver + lat_adm_device(:,ind_lat);
    lat_adm(ind_lat) = mean(cur_lat_adm_perserver);
    cur_lat_AR_perserver = cur_lat_AR_perserver + lat_adm_device(:,ind_lat);
    lat_allRed(ind_lat) = mean(cur_lat_AR_perserver);
end

figure;
semilogy(lat_allRed/1000,abs(util-opt)/opt,1:length(utildc),abs(utildc-opt)/opt,1:length(utiladm),abs(utiladm-opt)/opt)

figure;
subplot(1,2,1);
semilogy(1:length(util),abs(util-opt)/opt,1:length(utildc),abs(utildc-opt)/opt,1:length(utiladm),...
    abs(utiladm-opt)/opt)
subplot(1,2,2);
semilogy(lat_allRed/1000,abs(util-opt)/opt,1:length(utildc),abs(utildc-opt)/opt,1:length(utiladm),abs(utiladm-opt)/opt)
