%% create time series

n=2000; %time series length
x=zeros(n,2);
x(1,:)=rand(1,2);
alpha1=12.5/500;  %fast rotation per iteration
alpha2=5*(1+exp(1)/100)/500;  %slow rotation per iteration
delta=0.00075;  %switching probability per iteration
for i=1:n,
    x(i+1,:)=cylinder_freqswitch(x(i,:),delta,alpha1,alpha2);
end
figure(1)
plot(x(:,1),x(:,2)) %plot raw orbit

g=cos(2*pi*x(:,2));
figure(2)
plot(g,'.-')


%% load nice trajectory

load nicedata2000.mat
g=cos(2*pi*traj(:,2));

%% embed

lag=10;

dim=3;

figure(3)
if dim==2,
    X=[g(1:end-2*lag) g(1+lag:end-lag)];
    plot(X(:,1),X(:,2),'.-')
elseif dim==3,
    X=[g(1:end-2*lag) g(1+lag:end-lag) g(1+2*lag:end)];
    plot3(X(:,1),X(:,2),X(:,3),'.-')
    xlabel('g(x_t)');ylabel('g(x_{t+\pi/2})');zlabel('g(x_{t+\pi})')
  %  view(-60,-45)
elseif dim==4,
    X=[g(1:end-3*lag) g(1+lag:end-2*lag) g(1+2*lag:end-lag) g(1+3*lag:end)];
    scatter3(X(:,1),X(:,2),X(:,3),5,X(:,4),'filled')
    view(-60,-45)
end

%% calculate eigenvectors

n=size(X,1);
disp("making transfer operator")
tic
stepforward=1;  %previously with embedded time series data used stepforward=5 and epsilonfactor=4
epsfactor=1;  %when working with raw skew product output embedded as a cylinder use stepforward=?? and epsilonfactor=??

neighbourcount=25;

Ptuned=make_TO(X,stepforward,neighbourcount);
toc
disp("stochasticising")
tic
Ps=stochasticise(Ptuned); %row stochasticise.
Ps=sparse(Ps);
toc
disp("calculating eigenvalues")

numevals=10;
tic
[u,v]=eigs(Ps,numevals);
toc
   

%% create transfer operator

function Ptuned=make_TO(X,stepforward,neighboursize)

n=size(X,1);
D=pdist2(X,X);
nn=zeros(n,1);
for i=1:n,
    s=sort(D(:,i),'ascend');
    nn(i)=s(neighboursize);
end

S=pdist2(X(1:end-stepforward,:),X(1+stepforward:end,:));

Ptuned = exp((-S.^2)./(nn(1:end-stepforward)*nn(1+stepforward:end)'));

end




