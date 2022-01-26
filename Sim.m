clear;
S=10; % Number of end-users per operator
O=4; % Number of operators
N=S*O; % Number of users
Rtot= 0.5*N; % Total resource
x2=Rtot/N*ones(N,1); % Initial values of allocations for 2-layer
x3=Rtot/N*ones(S,O); % Initial values of allocations for 3-layer
x3u=Rtot/N*ones(S,O); % Initial used values of allocations for 3-layer
lambda=1; % Initial lambda
mu=lambda*ones(O,1); % Initial mu
y=Rtot/O*ones(O,1); % Initial y
yu=y; % Initial used y

% number of iterations and step sizes
k1=10;
k2=10;
g2=1/80;
g1=2;
ds_dm=1/80;

Niter=120;
k=Niter;

w=rand(N,1)+1; % initial weights
h = sqrt(1/2)*(randn(N,1)+1i*randn(N,1)); % initial channel
a = abs(h.^2);


T=(1/3)*(Niter); % duty cycle
p=lambda;

x2save=zeros(N,Niter);
ysave=zeros(O,Niter*k1);
ytemp=zeros(O,k1);
musave=zeros(O,Niter*k1*k2);
mutemp2=zeros(O,k1*k2);
mutemp=zeros(O,k2);
psave=zeros(1,Niter);
step2=1/160;

f2=zeros(1,Niter);
f2converge=zeros(1,Niter);
Res=zeros(1,Niter);
x2u=Rtot/N*ones(N,1);
x3u=Rtot/N*ones(S,O);
f3o=zeros(O,1);
f3oj=zeros(O,k1*k2*k);
f3c=zeros(O,1);
f3cj=zeros(O,k1*k2*k);

tol=10;
ptemp=0;
for i=1:Niter
    flag=0;
    if(mod(i,T)==0) % renew parameters each time slot
        w=rand(N,1)+1;
        h = sqrt(1/2)*(randn(N,1)+1i*randn(N,1));
        a = abs(h.^2);
        flag=1;
    end
    %     if(abs(p-ptemp)<tol)
    %         flag=1;
    %
    %     end
    w3=reshape(w,[S,O]);
    a3=reshape(a,[S,O]);
    x2=w./p - 1./a;
    x2=min(max(0,x2),Rtot);
    
    if sum(y)-Rtot > 0
        yu=y-((sum(y)-Rtot)/O).*ones(O,1);
        yu=min(max(0,yu),Rtot);
    else
        yu=y;
    end
    x2u=x2;
    if sum(x2)-Rtot > 0
        x2u=x2-((sum(x2)-Rtot)/N).*ones(N,1);
        x2u=min(max(0,x2u),Rtot);
        while sum(x2u)-Rtot > 0
        x2u=x2u-((sum(x2u)-Rtot)/N).*ones(N,1);
        x2u=min(max(0,x2u),Rtot);
        end
    else
        x2u=x2;
    end
    
    x2save(:,i)=x2;
    ptemp=p;
    p=p+step2*(sum(x2)-Rtot);
    p=max(0,p);
    psave(i)=p;
    
    Res(i)=sum(x2);
    Resu(i)=sum(x2u);
    f2(i)=sum(w.*log(1+a.*x2u));
    f2converge(i)=sum(w.*log(1+a.*x2));
    
    
    
    
    
    %% 3-Layer
    
    
    
    
    
    lambda=lambda+step2*(sum(y)-Rtot);
    lambda=max(0,lambda);
    
    for o=1:O
        for j=1:k1
            y(o)=y(o)+g1*(mu(o)-lambda);
            y(o)=max(0,y(o));
            ytemp(o,j)=y(o);
            
            %             if (j==2 && (flag==1 || i==1))
            %                 x3u(:,o)=x3(:,o);
            %             end
            for h=1:k2
                
                x3(:,o)=w3(:,o)./mu(o)-1./a3(:,o);
                x3(:,o)=min(max(0,x3(:,o)),Rtot);
                
                if sum(x3(:,o))-yu(o) > 0
                    x3u(:,o)=x3(:,o)-((sum(x3(:,o))-yu(o))/S).*ones(S,1);
                    x3u(:,o)=min(max(0,x3u(:,o)),Rtot);
                else
                    x3u(:,o)=x3(:,o);
                end
                
                mu(o)=mu(o)+g2*(sum(x3(:,o))-y(o));
                mu(o)=max(0,mu(o));
                f3o(o)=sum(w3(:,o).*log(1+a3(:,o).*x3u(:,o)));
                f3c(o)=sum(w3(:,o).*log(1+a3(:,o).*x3(:,o)));
                f3oj(o,(i-1)*k1*k2+(j-1)*k2+h)=f3o(o);
                f3cj(o,(i-1)*k1*k2+(j-1)*k2+h)=f3c(o);
                mutemp(o,h)=mu(o);
                
            end
            mutemp2(:,((j-1)*k2+1:(j*k2))) = mutemp;
        end
    end
    ysave(:,((i-1)*k1+1:(i*k1))) = ytemp;
    musave(:,((i-1)*k1*k2+1:(i*k1*k2))) = mutemp2;
    
end
f3=sum(f3oj);
f3converge=sum(f3cj);


scale = N/(O+k1*k2*ds_dm*S);

f2long1=zeros(round(scale*k1*k2),1);
f2convergelong1=zeros(round(scale*k1*k2),1);
for i=1:Niter
    f2long1(:,i)=f2(i)*ones(round(scale*k1*k2),1);
    f2convergelong1(:,i)=f2converge(i)*ones(round(scale*k1*k2),1);
end
f2long=reshape(f2long1,[1,round(scale*k1*k2)*Niter]);
f2convergelong=reshape(f2convergelong1,[1,round(scale*k1*k2)*Niter]);

f2longIter=zeros(k1*k2,1);
f2convergelongIter=zeros(k1*k2,1);
for i=1:Niter
    f2longIter(:,i)=f2(i)*ones(k1*k2,1);
    f2convergelongIter(:,i)=f2converge(i)*ones(k1*k2,1);
end
f2longIter=reshape(f2longIter,[1,k1*k2*Niter]);
f2convergelongIter=reshape(f2convergelongIter,[1,k1*k2*Niter]);

ScalingZeros=length(f2convergelong)/3-length(f3converge)/3;
f3convergeSC=[f3converge(1,1:k1*k2*(T-1)) f3converge(k1*k2*(T-1))*ones(1,round(scale*k1*k2)*(T-1)-k1*k2*(T-1)) f3converge(1,k1*k2*(T-1)+1:k1*k2*(2*T-1))...
    f3converge(k1*k2*(2*T-1))*ones(1,round(scale*k1*k2)*(T)-k1*k2*(T)) f3converge(1,k1*k2*(2*T-1)+1:k1*k2*(3*T-1)) ...
    f3converge(k1*k2*(3*T-1))*ones(1,round(scale*k1*k2)*(T)-k1*k2*(T)) ];

%fTh=[f2convergelong((T-1)*k1*k2)*ones(k1*k2*(T-1),1); f2convergelong((2*T-1)*k1*k2)*ones(k1*k2*(T),1); f2convergelong((3*T-1)*k1*k2)*ones(k1*k2*(T),1) ];

%plot((1:k1*k2*Niter)/(k1*k2),f3,(1:k1*k2*Niter)/(k1*k2),f2long,(1:k1*k2*Niter)/(k1*k2),f3converge,(1:k1*k2*Niter)/(k1*k2),f2convergelong,(1:k1*k2*Niter)/(k1*k2),fTh)
%%
thefg1=[f2convergelong(9438)*ones(1,9438) f2convergelong(19118)*ones(1,9680) f2convergelong(28798)*ones(1,9922)];
plot((1:round(scale*k1*k2)*Niter)/round(scale*k1*k2)*0.1,thefg1);
hold on
plot((1:length(f3convergeSC))/round(scale*k1*k2)*0.1,f3convergeSC,(1:round(scale*k1*k2)*Niter)/round(scale*k1*k2)*0.1,f2convergelong)
xlabel('Time');
ylabel('Total Utility');
hold off

figure;
thefg2=[f2convergelongIter(3900)*ones(1,3900) f2convergelongIter(7900)*ones(1,4000) f2convergelongIter(11900)*ones(1,4100)];
plot((1:k1*k2*Niter)/(k1*k2),thefg2);
hold on
plot((1:length(f3converge))/(k1*k2),f3converge,(1:k1*k2*Niter)/(k1*k2),f2convergelongIter)
xlabel('Top layer iterations $j$');
ylabel('Total Utility');
hold off
figure;
plot((1:k1*k2*Niter)/(k1*k2),f2convergelongIter,(1:k1*k2*Niter)/(k1*k2),f2longIter)

xlabel('Top layer iterations $j$');
ylabel('Total Utility');
%%
% figure;
% muopt=ones(1,100)*musave(4,100);
% errMu=abs(muopt-musave(4,(1:100)));
% secDer=((1+a.*x2).^2)./(w.*a.^2);
% rateMu=max(abs(1-g2*50),abs(1-g2*1));
% rateMu2=max(abs(1-0.1*g2*50),abs(1-0.1*g2*1));
% %rateMu=0.99;
% Hmu=10;
% t=(1:100000);
% TerrMu=Hmu.*rateMu.^t + ((1-rateMu.^t)./(1-rateMu)).*g2.*0.01;
% TerrMu2=Hmu.*rateMu2.^t + ((1-rateMu2.^t)./(1-rateMu2)).*0.1.*g2.*0.01;
% semilogy((1:100000),TerrMu,(1:100000),TerrMu2)
%%
% figure;
% yopt=[ysave(:,3898)*ones(O,3900) psave(7898)*ones(O,4000) psave(11900)*ones(1,4000)];
% err=abs(ysave-yopt);
% k=1:3900;
% k1=1:4000;
% k2=1:21;
% rate=(15*15-8)/(15*15+8);
% bound1=abs(psave(39)-psave(1)).*rate.^k./(1-rate.^k);
% bound2=abs(psave(79)-psave(39)).*rate.^k1./(1-rate.^k1);
% bound3=abs(psave(100)-psave(79)).*rate.^k2./(1-rate.^k2);
% therr=[bound1 bound2 bound3];
% semilogy((1:100)/10,err,(1:100)/10,therr)
% xlabel('$\|h\|$','interpreter','latex')

%%
% figure;
% popt=[psave(39)*ones(1,39) psave(79)*ones(1,40) psave(100)*ones(1,21)];
% err=abs(psave-popt);
% k=1:39;
% k1=1:40;
% k2=1:21;
% rate=(15*15-8)/(15*15+8);
% bound1=abs(psave(39)-psave(1)).*rate.^k./(1-rate.^k);
% bound2=abs(psave(79)-psave(39)).*rate.^k1./(1-rate.^k1);
% bound3=abs(psave(100)-psave(79)).*rate.^k2./(1-rate.^k2);
% therr=[bound1 bound2 bound3];
% semilogy((1:100)/10,err,(1:100)/10,therr)
% xlabel('$\|h\|$','interpreter','latex')