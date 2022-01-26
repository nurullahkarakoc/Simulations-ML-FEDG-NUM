%clear;
S=10; % Number of end-users per operator
O=4; % Number of operators
N=S*O; % Number of users
Rtot= 0.5*N; % Total resource
x2=Rtot/N*ones(N,1); % Initial values of allocations for 2-layer
x3=Rtot/N*ones(S,O); % Initial values of allocations for 3-layer
%x3u=Rtot/N*ones(S,O); % Initial used values of allocations for 3-layer
lambda=1.1; % Initial lambda
mu=lambda*ones(O,1); % Initial mu
y=Rtot/O*ones(O,1); % Initial y
yu=y; % Initial used y

% number of iterations and step sizes
k1=1;
k2=1;
g2=2/52.5;
g1=4;
ds_dm=1/80;

Niter=100;
k=Niter;

w=rand(N,1)+1; % initial weights
%h = sqrt(1/2)*(randn(N,1)+1i*randn(N,1)); % initial channel
%a = abs(h.^2);
a=ones(N,1);
alphabar=5;
betabar=0.25;

T=(1)*(Niter); % duty cycle
p=lambda;

x2save=zeros(N,Niter);
ysave=zeros(O,Niter*k1);
ytemp=zeros(O,k1);
musave=zeros(O,Niter*k1*k2);
mutemp2=zeros(O,k1*k2);
mutemp=zeros(O,k2);
psave=zeros(1,Niter);
lambdasave=zeros(1,Niter);
step2=2/210;

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


for i=1:Niter-1
    flag=0;
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
    lambdasave(i)=lambda;
    
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
            mutemp2(o,((j-1)*k2+1:(j*k2))) = mutemp(o,:);
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

ScalingZeros=length(f2convergelong)/4-length(f3converge)/4;
%f3convergeSC=[f3converge(1,1:k1*k2*(T-1)) f3converge(k1*k2*(T-1))*ones(1,round(scale*k1*k2)*(T-1)-k1*k2*(T-1)) f3converge(1,k1*k2*(T-1)+1:k1*k2*(2*T-1))...
%    f3converge(k1*k2*(2*T-1))*ones(1,round(scale*k1*k2)*(T)-k1*k2*(T)) f3converge(1,k1*k2*(2*T-1)+1:k1*k2*(3*T-1)) ...
%    f3converge(k1*k2*(3*T-1))*ones(1,round(scale*k1*k2)*(T)-k1*k2*(T)) f3converge(1,k1*k2*(3*T-1)+1:k1*k2*(4*T-1)) ...
%    f3converge(k1*k2*(4*T-1))*ones(1,round(scale*k1*k2)*(T)-k1*k2*(T))];

%fTh=[f2convergelong((T-1)*k1*k2)*ones(k1*k2*(T-1),1); f2convergelong((2*T-1)*k1*k2)*ones(k1*k2*(T),1); f2convergelong((3*T-1)*k1*k2)*ones(k1*k2*(T),1) ];

%plot((1:k1*k2*Niter)/(k1*k2),f3,(1:k1*k2*Niter)/(k1*k2),f2long,(1:k1*k2*Niter)/(k1*k2),f3converge,(1:k1*k2*Niter)/(k1*k2),f2convergelong,(1:k1*k2*Niter)/(k1*k2),fTh)
%plot((1:length(f3convergeSC))/round(scale*k1*k2)*0.1,f3convergeSC,(1:round(scale*k1*k2)*Niter)/round(scale*k1*k2)*0.1,f2convergelong)
%xlabel('Time');
%ylabel('Total Utility');
figure;
plot((1:length(f3converge))/(k1*k2),f3converge,(1:k1*k2*Niter)/(k1*k2),f2convergelongIter)
figure;
plot((1:k1*k2*Niter)/(k1*k2),f2convergelongIter,(1:k1*k2*Niter)/(k1*k2),f2longIter)

xlabel('Top layer iterations $j$');
ylabel('Total Utility');
%%
iter=1:99;
ratelambda=(5-0.25)/(5+0.25);
boundlambda=abs(1.1-psave(99))*(ratelambda.^iter);

err=abs(sum(x2save)-sum(ysave(:,1:1:100)));
eps=max(err(1:99));
errBoundlambda=abs(1-psave(99))*(ratelambda.^iter)+step2*eps*(1-ratelambda.^iter)./(1-ratelambda);
%%
figure;
plot(iter,abs(psave(1:99)-psave(99))/psave(99),iter,boundlambda);
xlabel('Iterations $j$');
ylabel('|$\lambda - \lambda^*|/|$\lambda^*$|');
hold on
%plot(iter,abs(lambdasave(1:99)-lambdasave(99))/lambdasave(99),iter,errBoundlambda);
hold on

%plot(iter,abs(lambdasave10(1:99)-lambdasave10(99))/lambdasave10(99),iter,errBoundlambda10);





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
% 
