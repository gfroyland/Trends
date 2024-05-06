function P=stochasticise(Pbad),

rowsum=sum(Pbad');
P=diag(sparse(1./rowsum))*Pbad;
