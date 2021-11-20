%% Model Parameters
alpha = 0.2;
beta = 0.5;
gamma = 0.9;
sigma = 0.1;
delta = 0.9;
%% Compute / Discretize Shock Distribution
nshocks =3;
[e,w] = qnwlogn(nshocks,-sigma^2/2,sigma^2);
%% Pack Model Structure
clear model
model.func ='myfunc';
model.discount = delta;
model.e = e;
model.w = w;
model.params = {alpha beta gamma};
%% Define Approximation Space
n = 10;
smin = 5;
smax = 10;
fspace = fundefn('cheb',n,smin,smax);
snodes = funnode(fspace);
estar = 1;
xstar = ((1-delta*gamma)/(delta*beta))^(1/beta-1);
sstar = gamma*xstar + xstar^beta;
pstar = (sstar-xstar).^(-alpha);
%% Compute LQ-Approximation
[vlq,xlq] = lqapprox(model,snodes,sstar,xstar,pstar);
%% Compute Bellman Equation
v = vlq; x = xlq;
Phi = funbas(fspace);
c = Phi\v;
maxit = 50; tol = 5e-8;
maxiti = 500; toli = sqrt(eps);
% solve via Newton (policy) iteration
for it = maxiti
    cold = c
    [v,x,vjac] = myvmax(snodes,x,c,e,maxit,tol,fspace,w,alpha,beta,gamma,delta)
    c = c-(Phi-vjace)\(Phi*c-v);
    if norm(c-cold)<toli, break, end
end
%% Diagnoses
% Diagnose if VF approximant solves BE:
nplot = 100;
splot = nodeunif(nplot,smin,smax);
% rough guess for actions at evaluation points:
x = funeval(Phi\x,fspace,splot);
% Values and actions at evaluation points:
[v,x]=mymax(splot,x,c,e,maxit,tol,fspace,w,alpha,beta,gamma,delta);
resid = v - funeval(c,fspace,splot);
plot(splot,resid)
% diagnose poor resid by bottom of page 236
for k=1:length(e)
    g = myfunc('g',snodes,x,e(k),alpha,beta,gamma);
    if any(g<smin), disp('Warning:_extrapolating_beyond_smin'), end;
    if any(g>smax), disp('Warning:_extrapolating_beyond_smax'), end;
end

% residual small and no warnings means that the VF approximant solves BE
