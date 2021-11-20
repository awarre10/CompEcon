%% Parameters
dbar = 1.0;
gamma = 0.5;
beta = 0.4;
sigma = 0.1;
delta = 0.9;
%% Discretize Shock Distribution
nshocks = 3;
[e,w] = qnwnorm(nshocks,0,sigma^2);
%% Basis Functions and Colocation Nodes
n = 10;
dmin = dbar+min(e)/(1-gamma);
dmax = dbar+max(e)/(1-gamma);
fspace = fundefn('cheb',n,dmin,dmax);
dnode = funnode(fspace);
%% Solving Decision Model
B = diag(dnode.^(-beta))*funbas(fspace,dnode);
y = 0;
for k=1:K
    dnext = dbar + gamma*(dnode-dbar) + e(k);
    B = B - delta*w(k)*diag(dnext.^(-beta))*funbas(fspace,dnext);
    y = y + delta*w(k)*dnext.^(1-beta);
end
c = B\y;
%% Response Function and Approximation Residuals
d = nodeunif(10*n,dmin,dmax);
p = funeval(c,fspace,d);
Eh=0;
for k=1:m
    dnext = dbar + gamma*(d-dbar) + e(k);
    h = diag(dnext.^(-beta))*(funeval(c,fspace,dnext)+dnext);
    Eh = Eh + delta*w(k)*h;
end
resid = d.^(-beta).*p-Eh;
plot(p,resid)