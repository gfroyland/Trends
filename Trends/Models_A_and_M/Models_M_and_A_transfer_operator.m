%% mean drift

step=0.1
t=0:step:100;
g=cos(t)+t/20;   %linear

%g=cos(t)+(t-50).^2/400;  %symmetric quadratic
%g=cos(t)+(t-40).^2/400;  %semi-symmetric quadratic
%g=cos(t)+(t-30).^2/400;  %nonsymmetric quadratic

%g(1:500)=cos(t(1:500))+(t(1:500)-50).^2/400; g(501:1001)=cos(t(501:end))-(t(501:end)-50).^2/400; %flipped quadratic
%g=g-t/20;

%g(1:500)=cos(t(1:500))+t(1:500)/20; g(501:1001)=cos(t(501:end))+t(501:end)/10-t(500)/20;   %piecewise linear 
%g(1:500)=cos(t(1:500))+t(1:500)/10; g(501:1001)=cos(t(501:end))+t(501:end)/20+t(500)/20; 

figure;
plot(g)

% embed in 3d
dim=3;
lag=round(pi/2/step);
X=[g(1:end-2*lag); g(1+lag:end-lag); g(1+2*lag:end)]';
figure
plot3(X(:,1),X(:,2),X(:,3))

%% amplitude drift

step=0.1;
t=0:step:100;
%g=(1+t/20).*cos(t);   dim=2;  
g=(1+5*sin(t/100*pi)).*cos(t);    dim=3

lag=pi/2/step;
figure;
plot(g,'.')
% embed in 
figure
if dim==2,
    X=[g(1:end-2*lag); g(1+lag:end-lag)]';
    plot(X(:,1),X(:,2),'.')
elseif dim==3,
    X=[g(1:end-2*lag); g(1+lag:end-lag); g(1+2*lag:end)]';
    plot3(X(:,1),X(:,2),X(:,3))
end


%% calc evecs and plot

stepforward=1;

epsilonfactor=1;  

neighbourcount=7;
numevals=15;
n=size(X,1);
[Ptuned]=make_TO(X,stepforward,neighbourcount);


Ps=stochasticise(Ptuned); %row stochasticise.
Ps=sparse(Ps);

%Ps=P;
[u,v]=eigs(Ps,numevals);
[ul,vl]=eigs(Ps',numevals);

%number of eigenfunction to be plotted below.
evnum=2;  %need 6 for linear amplitude drift, 8 for quadratic amplitude drift, 2 for all mean drifts
if evnum==1,
    ul(:,evnum)=ul(:,evnum)*sign(ul(end,evnum))*sqrt(size(ul,1));
    figure
    subplot(1,3,1);scatter(t(1+lag:end-lag-stepforward),X(1:end-stepforward,1),3,ul(:,evnum),'filled'); xlabel('time');ylabel('g(x_t)');colorbar;colormap jet
    subplot(1,3,2);plot(t(1+lag:end-lag-stepforward),ul(:,evnum)); xlabel('time'); ylabel('eigenvector value')
    subplot(1,3,3);scatter3(X(1+lag:end-lag-stepforward,1),X(1+lag:end-lag-stepforward,2),X(1+lag:end-lag-stepforward,3),5,ul(1+lag:end-lag,evnum),'filled');xlabel('g(x_t)');ylabel('g(x_{t+\pi/2})');zlabel('g(x_{t+\pi})');colorbar;colormap jet
else
    u(:,evnum)=u(:,evnum)*sign(u(end,evnum))
    u(:,evnum)=u(:,evnum)/norm(u(:,evnum))*sqrt(size(u,1));
    figure
 subplot(1,3,1);scatter(1:size(X,1)-stepforward,X(1:end-stepforward,1),3,u(:,evnum),'filled'); xlabel('time');ylabel('g(x_t)');colorbar;colormap jet
    subplot(1,3,2);plot(1:size(u,1),u(:,evnum)); xlabel('time'); ylabel('eigenvector value')
     if dim==2,
        subplot(1,3,3);scatter(X(1+lag:end-lag-1,1),X(1+lag:end-lag-1,2),5,u(1+lag:end-lag,evnum),'filled');xlabel('g(x_t)');ylabel('g(x_{t+\pi/2})');zlabel('g(x_{t+\pi})');colorbar;colormap jet
    elseif dim==3,
        subplot(1,3,3);scatter3(X(1+lag:end-lag-stepforward,1),X(1+lag:end-lag-stepforward,2),X(1+lag:end-lag-stepforward,3),5,u(1+lag:end-lag,evnum),'filled');xlabel('g(x_t)');ylabel('g(x_{t+\pi/2})');zlabel('g(x_{t+\pi})');colorbar;colormap jet
    end
 end    


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




