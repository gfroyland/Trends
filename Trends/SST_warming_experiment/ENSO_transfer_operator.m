lag=12;

load IP_SSTA_ERSSTv4-large.mat 

x=ssta;

X=[x(:,1:end-2*lag); x(:,1+lag:end-lag)];

dim=2;
X=[];
for d=1:dim,
    X=[X; x(:,1+(d-1)*lag:end-(dim-d)*lag)];
end
X=X';


%% make TO

L2=1;  %1 if L2 norm, 0 if wasserstein distance
n=size(X,1);
disp("making transfer operator")
tic
stepforward=6; 
neighboursize=5;

Ptuned = make_TO(X,stepforward,neighboursize);
toc

%% calc evecs and plot
disp("stochasticising")
tic
Ps=stochasticise(Ptuned); %row stochasticise.
Ps=sparse(Ps);
toc
disp("calculating eigenvalues")

numevals=15;
evnum=2;
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



