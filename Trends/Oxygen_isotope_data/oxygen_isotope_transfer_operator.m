%% load data

N=3000;
load d18O.mat
g=Obs.d18O(1:end-1);


%% embed

lag=15;

dim=5;
X=[];
for d=1:dim,
    X=[X g(1+(d-1)*lag:end-(dim-d)*lag)];
end

%% calculate eigenvectorss 

n=size(X,1);
disp("making transfer operator")
tic
stepforward=7;  %previously with embedded time series data used stepforward=5 and epsilonfactor=4
neighboursize=7;
%epsfactor=1;  %when working with raw skew product output embedded as a cylinder use stepforward=?? and epsilonfactor=??
Ptuned=make_TO(X,stepforward,neighboursize);
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

