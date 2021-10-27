% PS #7: Andres Warren %

% Exercise 1 %

% 5 degree Chebychev %
alpha = 1;
fspace5 = fundefn('cheb',5,-1,1);
c5 = funfitf(fspace5,'f7',alpha);
x5 = nodeunif(1001, -1,1);
yact5 = f7(x5, alpha);
yapp5 = funeval(c5, fspace5, x5);
plot(x5, yact5-yapp5)
title('Plot 1')

% 50 degree Chebychev %
alpha = 1;
fspace50 = fundefn('cheb',50,-1,1);
c50 = funfitf(fspace50,'f7',alpha);
x50 = nodeunif(1001, -1,1);
yact50 = f7(x55, alpha);
yapp50 = funeval(c55, fspace50, x50);
plot(x50, yact50-yapp50)
title('Plot 2')

% 5 degree linear spline %
alpha = 1;
fspace5lin = fundefn('spli',5,-1,1,1);
c5lin = funfitf(fspace5lin,'f7',alpha);
x5lin = nodeunif(1001, -1,1);
yact5lin = f7(x5lin, alpha);
yapp5lin = funeval(c5lin, fspace5lin, x5lin);
plot(x5lin, yact5lin-yapp5lin)
title('Plot 3')

% 50 degree linear spline %
alpha = 1;
fspace50lin = fundefn('spli',50,-1,1,1);
c50lin = funfitf(fspace50lin,'f7',alpha);
x50lin = nodeunif(1001, -1,1);
yact50lin = f7(x50lin, alpha);
yapp50lin = funeval(c50lin, fspace50lin, x50lin);
plot(x50lin, yact50lin-yapp50lin)
title('Plot 4')

% 5 degree cubic spline %
alpha = 1;
fspace5cub = fundefn('spli',5,-1,1);
c5cub = funfitf(fspace5cub,'f7',alpha);
x5cub = nodeunif(1001, -1,1);
yact5cub = f7(x5cub, alpha);
yapp5cub = funeval(c5cub, fspace5cub, x5cub);
plot(x5cub, yact5cub-yapp5cub)
title('Plot 5')

% 50 degree cubic spline %
alpha = 1;
fspace50cub = fundefn('spli',50,-1,1);
c50cub = funfitf(fspace50cub,'f7',alpha);
x50cub = nodeunif(1001, -1,1);
yact50cub = f7(x50cub, alpha);
yapp50cub = funeval(c50cub, fspace50cub, x50cub);
plot(x50cub, yact50cub-yapp50cub)
title('Plot 6')

% Exercise 2 %

alpha = 1; eta = 1.5;
n = 25; a = 0.3; b = 3;
fspace = fundefn('cheb', n, a,b);
p = funnode(fspace);
Phi = funbas(fspace, p);
c = Phi\sqrt(p);
c = broyden('cournot_resid', c, p, alpha, eta, Phi);

% supply curves %
plot(funeval(c, fspace, p),p)
hold on
plot(3*funeval(c, fspace, p),p)
plot(5*funeval(c, fspace, p),p)
plot(10*funeval(c, fspace, p),p)
plot(20*funeval(c, fspace, p),p)

% demand curve %
plot(p.^(-eta),p)

% Excercise 3 %

% From the book %
r = 0.1;
k = 0.5;
eta = 5;
s0 = 1;

T = 1;
n = 15;
tnodes = chebnode(n-1,0,T);
fspace = fundefn('cheb',n,0,T);

c = zeros(n,2); c(1,:) = 1;
c = broyden('resid',c(:),tnodes,T,n,fspace,r,k,eta,s0);

nplot 501;
t = nodeunif(nplot,0,T);
c = reshape(c,n,2);
x = funeval(c,fspace,t);
d = funeval(c,fspace,t,1);
R = d - [r*x(:,1)+k -x(:,1).^(-eta)];


% Define model
 A=[-1 -.5;0 -.5];
 model.func='pbvp01';
 model.tb=[0;1];
 model.params={A};

% Define approximant
 n=5;
 a=0;
 b=2;
 fspace=fundefn('cheb',n-1,a,b);
 tnodes=funnode(fspace);
 fspace=fundefn('cheb',n,a,b);

% Initial conditions
 c=zeros(fspace.n,2);

% Evaluation points for plots
 tvals=linspace(a,b,201)';

% Call solver
 optset('broyden','defaults')
 [c,x,r]=bvpsolve(model,fspace,tnodes,c,tvals);

% Produce plots
 close all
 figure(1)
 C=exp(0.5)*(1-exp(-1));
 plot(tvals,exp(-tvals)-x(:,1),'-', ...
     tvals,C*exp(-tvals/2)+exp(-tvals)-x(:,2),'--');
 title('BVP Example: Approximation Errors')
 xlabel('t')
 ylabel('x-\phi(t)c')
 legend('x_1','x_2')

 figure(2)
 plot(tvals,r(:,1),'-',tvals,r(:,2),'--');
 title('BVP Example: Residual Functions')
 xlabel('t')
 ylabel('r')
 legend('x_1','x_2')
 
