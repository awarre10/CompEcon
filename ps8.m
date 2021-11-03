% Problem Set #8 %

% Exercise 1 %

% C:
T = 0.75;
sigma = 0.15;
r = 0.06;
strike = 1.75;
p0 = 1.50;
N = 200;
tau = T/N;
delta = exp(-r*tau);
u = exp(sigma*sqrt(tau));
q = 0.5+sqrt(tau)*(r-(sigma^2)/2)/(2*sigma);
price = p0*(u.^(-N:N))';
n = length(price);
f = [zeros(n,1) strike-price];
P = zeros(2,n,n);
for i=1:n
P(1,i,min(i+1,n)) = q;
P(1,i,max(i-1,1)) = 1-q;
end
model.reward = f;
model.transprob = P;
model.horizon = N;
model.discount = delta;
[v,x] = ddpsolve(model);
% D:
plot(price, v(:,1)); axis([0 strike*2 -inf inf])

% Exercise 2 %

% C:
sigma = 0.2;                                          
T     = 0.5;                                          
K     = 1.0;                                          
p0    = 1.0;                                          
r     = 0.1;                                          
nt = 30;
deltat = T/(nt-1);
beta = exp(-r*deltat);
m = 10;
[e,w] = qnwnorm(m,deltat*(r-sigma^2/2),deltat*sigma^2);
x = [0;1];
nx = length(x);
clear model;
model.func = 'american_option';
model.horizon = nt;
model.discount = beta;
model.e = e;
model.w = w;
model.actions = x;
model.discretestates = 2;
model.params = {K};
d = 6;
ns = 500;
s0 = log(p0);
smin = s0-d*sigma*sqrt(T);
smax = s0+d*sigma*sqrt(T);
fspace = fundefn('lin', ns, smin, smax,[],[0;1]);
scoord = funnode(fspace);
s = gridmake(scoord);
p = exp(s(:,1));
v = zeros(2*ns,1);
optset('dpsolve','nres',5);
[c,s,v,x] = dpsolve(model, fspace, s,v);

% D:
%figure(1);
p = exp(scoord{1});
v = reshape(v,ns,2,nt+1);
v = squeeze(v(:,1,:));
plot(p,v(:,1),p,max(K-p,0));
title('American Put Option Value');
xlabel('Asset Price'); ylabel('Value');

% E:
%figure(2);
x = reshape(x,ns,2,nt);
x = squeeze(x(:,1,:)); 
plot(linspace(T,0,nt),p(sum(x==1)))
title('American Put Option Optimal Exercise Boundary'); 
xlabel('Time-to-Maturity'); ylabel('Asset Price');

% Exercise 3

% C:
delta = 0.9;alpha = 0.2;beta = 0.5;gamma = 0.9;sigma = 0.1;
m = 3;
[e,w] = qnwlogn(m,0,sigma^2);

n = 10;
smin = 5;
smax = 10;
fspace = fundefn('cheb',n,smin,smax);
snodes = funnode(fspace);

model.func = 'optimal_growth';
model.discount = delta;
model.e = e;
model.w = w;
model.params = {alpha beta gamma};

estar = 1;
xstar = ((1-delta*gamma)/(delta*beta))^(1/(beta-1));
sstar = gamma*xstar + xstar^beta;

[vlq, xlq] = lqapprox(model, snodes, sstar, xstar, estar);

[c,s,v,x,resid] = dpsolve(model,fspace,snodes,vlq,xlq);

% D: the utility derived from consuming a unit of the good today must
% equal the discounted expected utility derived from investing it and consuming its yield tomorrow.

% E: themarginal product of capital must equal the capital depreciation rate plus the discount rate.

% F:
plot(s, x./s, snodes,xlq./snodes)

% G:
nyrs = 20;
nrep = 2000;
sinit = 5*ones(nrep,1);
[spath, xpath] = dpsimul(model, sinit, nyrs,s,x);

nsmooth = 5;
nbin = 80;
[ss,pi,xx] = dpstst(model, nsmooth, nbin, s,x);

plot(ss,pi)
title("wealth distribution for steady state")

% H:
plot(linspace(1,nyrs,nyrs+1),mean(spath))
title("Expected wealth over time")

% Exercise 4 %

% E:
alpha = [0.8 0.2];
beta = [0.8 0.2;0.13 0.5];
gamma = [-0.8 0];
omega = [0.3 1];
starget = [0 1];
sigma = 0.04*eye(2);
delta = 0.9;

m = [5 5];
mu = [0 0];
[e,w] = qnwnorm(m,mu,sigma);

model.func = 'monetary_policy';
model.discount = delta;
model.e = e;
model.w = w;
model.params = {alpha beta gamma omega starget};

n = [20 20];
nn = prod(n);
smin = [-15 -10];
smax = [15 10];
fspace = fundefn('spli',n, smin, smax);
scoord = funnode(fspace);
snodes = gridmake(scoord);

estar = mu;
sstar = starget;
xstar = (sstar(1)-alpha(1)-beta(1,:)*sstar')/gamma(1);
pstar = [0 0];

[vlq,xlq,plq,ss,xx,pp] = lqapprox(model,snodes,sstar,xstar,pstar);

optset('dpsolve','nres',4);
[c,s,v,x] = dpsolve(model, fspace, snodes,vlq, xlq);

% F:
figure(1)
hh = surf(s{1},s{2},x');
title('Optimal Monetary Policy');
xlabel('GDP Gap');
ylabel('Inflation Rate');
zlabel('Norminal Interest Rate');

% G:
figure(2)
contour(s{1},s{2},x',25);
title('Optimal Monetary Policy');
xlabel('GDP Gap');
ylabel('Inflation Rate');
zlabel('Norminal Interest Rate');
rotate3d

% H:
nyrs = 20;
nrep = 5000;
sinit = smax(ones(nrep,1),:)-1;
[spath,xpath] = dpsimul(model, sinit, nyrs, s,x);
s1path = squeeze(spath(:,1,:));
s2path = squeeze(spath(:,2,:));

figure(3);
plot(0:nyrs,mean(s1path));
title('Expected State Path');
xlabel('Year');
ylabel('GDP Gap');

% I:
figure(4);
plot(0:nyrs,mean(s2path));
title('Expected State Path');
xlabel('Year');
ylabel('Inflation');

% Exercise 6 %

% Rational Expectations Model: f(s,x(s),E[h(g(s,x(s),ϵ),x(g(s,x(s),ϵ)))])=0 %

% State Transition: st+1g(st,xt,ϵt+1) %

% ENTER MODEL PARAMETERS
  dbar  = 1.0;                                          % long-run mean dividend
  gamma = 0.5;                                          % dividend autoregression coefficient
  beta  = 0.5;                                          % coefficient of risk aversion
  sigma = 0.1;                                          % dividend volatility
  delta = 0.9;                                          % discount factor

% COMPUTE SHOCK DISTRIBUTION
  m = 3;                                                % number of production shocks
  [e,w] = qnwnorm(m,0,sigma^2);                         % normal nodes and weights
  
% DEFINE APPROXIMATION SPACE
  n      = 10;                                          % degree of approximation
  dmin   = dbar+min(e)/(1-gamma);                       % minimum production
  dmax   = dbar+max(e)/(1-gamma);                       % maximum production
  fspace = fundefn('cheb',n,dmin,dmax);                 % function space
  dnode  = funnode(fspace);                             % state collocaton nodes

% PACK MODEL STRUCTURE
  clear model
  model.func = 'mfrem01';                               % model functions
  model.discount = delta;                               % discount factor
  model.e = e;                                          % shocks
  model.w = w;                                          % probabilities
  model.params = {delta dbar gamma beta};               % other parameters
  
% INITIALIZE RESPONSE
  xinit = dbar/(1-delta)+(dnode-dbar)/(1-gamma*delta); 
  
% SOLVE RATIONAL EXPECTATIONS EQULIBRIUM
  optset('remsolve','nres',10);
  optset('remsolve','showiters',0);
tic
  [c,d,p,x,f,resid] = remsolve(model,fspace,dnode,xinit); 
 toc
% PLOT EQUILIBRIUM PRICE
  figure(1); 
  plot(d,p);
  title('Equilibrium Pricing Function')
  xlabel('Dividend Level'); 
  ylabel('Asset Price');
  
% PLOT RESIDUAL
  figure(2); 
  plot(d,resid);
  title('Approximation Residual')
  xlabel('Dividend Level'); 
  ylabel('Residual');
  
% PLOT ARBITRAGE PROFIT
  figure(3); 
  plot(d,f);
  title('Arbitrage Profit Function')
  xlabel('Dividend Level'); 
  ylabel('Arbitrage Profit');
  
% SOLVE RATIONAL EXPECTATIONS EQULIBRIUM - DIRECT METHOD
  LHS = diag(dnode.^(-beta))*funbas(fspace,dnode);
  RHS = 0;
  for k=1:m
    dnext = dbar + gamma*(dnode-dbar) + e(k);
    LHS   = LHS - delta*w(k)*diag(dnext.^(-beta))*funbas(fspace,dnext);
    RHS   = RHS + delta*w(k)*dnext.^(1-beta);
  end
  c = LHS\RHS;

% COMPUTE RESIDUAL
  d = nodeunif(10*n,dmin,dmax);
  p = funeval(c,fspace,d);
  Eh=0;
  for k=1:m
    dnext = dbar + gamma*(d-dbar) + e(k);
    h     = diag(dnext.^(-beta))*(funeval(c,fspace,dnext)+dnext);
    Eh    = Eh + delta*w(k)*h;
  end
  resid = d.^(-beta).*funeval(c,fspace,d)-Eh;

% PLOT EQUILIBRIUM PRICE
  figure(4); 
  plot(d,p);
  title('Equilibrium Pricing Function')
  xlabel('Dividend Level'); 
  ylabel('Asset Price');
  
% PLOT RESIDUAL
  figure(5); 
  plot(d,resid);
  title('Approximation Residual')
  xlabel('Dividend Level'); 
  ylabel('Residual');
