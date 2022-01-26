%clear;
%S=20; % Number of end-users per operator
%O=10; % Number of operators
%N=S*O; % Number of users
N=O;
O=L;
Rtot= 0.5*N; % Total resource
x2=Rtot/N*ones(N,1); % Initial values of allocations for 2-layer
x3=Rtot/N*ones(N,1); % Initial values of allocations for 3-layer
x3dc=Rtot/N*ones(N,1); % Initial used values of allocations for 3-layer

lambda=1; % Initial lambda
mu=lambda*ones(L,1); % Initial mu
mudc=mu;
lambdadc=mu;
y=Rtot/L*ones(L,1); % Initial y
ydc=y; % Initial used y

% W=erdos_renyi(O,0.35);
% idx = W ~=0;
% c = sum(idx,1);
% max_Neighbor = max(c)-1;
% min_Neighbor = min(c)-1;
%     while min_Neighbor == 0
%         W=erdos_renyi(O,0.25);
%         idx = W ~=0;
%         c = sum(idx,1);
%         max_Neighbor = max(c)-1;
%         min_Neighbor = min(c)-1;
%     end

% number of iterations and step sizes
k1=3;
k2=3;
g2=2/105;
g1=10;
step2=2/1000;
ds_dm=1/5;
phi=1;

Niter=12000;
k=Niter;

w=rand(N,1)+1; % initial weights
h = sqrt(1/2)*(randn(N,1)+1i*randn(N,1)); % initial channel
a = abs(h.^2);

T=(1)*(Niter); % duty cycle
p=lambda;

x2save=zeros(N,Niter);
ysave=zeros(O,Niter*k1);
ydcsave=zeros(O,Niter*k1);
ytemp=zeros(O,k1);
ydctemp=zeros(O,k1);
musave=zeros(O,Niter*k1*k2);
mudcsave=zeros(O,Niter*k1*k2);
mutemp2=zeros(O,k1*k2);
mudctemp2=zeros(O,k1*k2);
mutemp=zeros(O,k2);
mudctemp=zeros(O,k2);
psave=zeros(1,Niter);
lambdasave=zeros(1,Niter);
lambdadcsave=zeros(O,Niter);


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
 %   w3=reshape(w,[S,O]);
 %   a3=reshape(a,[S,O]);
  %  x2=w./p - 1./a;
  %  x2=min(max(0,x2),Rtot);
    
%     if sum(y)-Rtot > 0
%         yu=y-((sum(y)-Rtot)/O).*ones(O,1);
%         yu=min(max(0,yu),Rtot);
%     else
%         yu=y;
%     end
  %  x2u=x2;
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
   % Resu(i)=sum(x2u);
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
    
    %% 3-Layer DC
    
    
    
    
    % Write consensus here
    lambdadc=lambdadc+(3*step2/i)*(O*ydc-Rtot);
    lambdadc=max(0,lambdadc);
    
    
    lambdadc=(W^phi)*lambdadc;
    lambdadcsave(:,i)=lambdadc;
    for o=1:O
        for j=1:k1
            ydc(o)=ydc(o)+g1*(mudc(o)-lambdadc(o));
            ydc(o)=max(0,ydc(o));
            ydctemp(o,j)=ydc(o);
            
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
                
                mudc(o)=mudc(o)+g2*(sum(x3dc(:,o))-ydc(o));
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
    ydcsave(:,((i-1)*k1+1:(i*k1))) = ydctemp;
    mudcsave(:,((i-1)*k1*k2+1:(i*k1*k2))) = mudctemp2;
    
end

f3=sum(f3oj);
f3converge=sum(f3cj);

f2longIter=zeros(k1*k2,1);
f2convergelongIter=zeros(k1*k2,1);
for i=1:Niter
    f2longIter(:,i)=f2(i)*ones(k1*k2,1);
    f2convergelongIter(:,i)=f2converge(i)*ones(k1*k2,1);
end
f2longIter=reshape(f2longIter,[1,k1*k2*Niter]);
f2convergelongIter=reshape(f2convergelongIter,[1,k1*k2*Niter]);
%%
plot((1:length(f3converge))/(k1*k2),f3converge,(1:k1*k2*Niter)/(k1*k2),f2convergelongIter,(1:length(f3))/(k1*k2),f3);


idx = W ~=0;
c = sum(idx,1);
max_Neighbor = max(c)-1;
scale1= N/(O+k1*k2*S*ds_dm);
scale2= N/(max_Neighbor+k1*k2*S*ds_dm);

figure;
plot((1:length(f3converge))/(k1*k2*scale1),f3converge,(1:k1*k2*Niter)/(k1*k2),f2convergelongIter,(1:length(f3))/(k1*k2*scale2),f3);

figure;
plot((1:Niter),lambdadcsave)

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