%% 11.1

% Define parameters
T=30;
kappa=.1;
alpha=.05;
sigma=0.1;

% Evaluate for values of the short rate between 0 and 2
%   (2 is as close to infinity as one needs to get for interest rates
%      unless the economy is in hyperinflation.)
n=20;
fspace=fundefn('cheb',n,0,2);

if 1
  % Use finsolve
  clear model
  model.func='mffin01';
  model.T=T;
  model.american=0;
  model.params={kappa,alpha,sigma};
  s=funnode(fspace);
  c=finsolve(model,fspace,'lines',s,1);
else
  % use the specialized function cirbond
  c=cirbond(fspace,T,kappa,alpha,sigma);
end

% Compute exact solution and plot solution and errors
sigma2=sigma*sigma;
gam=sqrt(kappa.^2+2*sigma2);
egam=exp(gam*T)-1;
denom=(gam+kappa).*egam+2*gam;
B=2*egam./denom;
A=(2*gam*exp((kappa+gam)*T/2)./denom).^(2*kappa*alpha./sigma2);

x=linspace(0,0.25,101)';
V=A*exp(-B*x);
Vhat=funeval(c(:,end),fspace,x);

figure(1)
plot(x,Vhat); 
title('Bond Price')
xlabel('r')
ylabel('V(r)')

figure(1)
plot(x,V-Vhat)
title('Approximation Errors')
xlabel('r')
ylabel('Error')

%% 11.2

% Define parameters
r      = 0.05;               % risk free interest rate
deltaS = 0;                  % dividend rate 
sigma  = 0.2;                % volatility  
K      = 1;                  % exercise price 
put    = 0;                  % put indicator
T      = 1;                  % time-to-maturity

% Create model variable
clear model
model.func='mffin02';
model.T=T;
model.params={r,deltaS,sigma,K,put};
model.american=0;

% Define approximation space
n=51; 
fspace=fundefn('lin',n,0,2*K);    
s=funnode(fspace);

% Call solution algorithm
c=finsolve(model,fspace,'lines',s);

% Compute exact solution and plot solution and errors
S=linspace(0,2*K,501)';
premium=bs(sigma,S,K,r,deltaS,T,put);

figure(2)
plot(S,funeval(c,fspace,S));
title('Call Option Premium')
xlabel('S')
ylabel('Premium')

figure(2)
plot(S,(premium-funeval(c,fspace,S)))
title('Approximation Errors')
xlabel('S')
ylabel('Error')

% Compute performance comparisons

stype={'lines' 'explicit' 'implicit' 'CN' 'stiff'};
N=[1 75 75 75 75];
m=length(stype);
C=zeros(n,m);
tl=zeros(1,m);

fspace=fundefn('lin',n,0,2*K);    
s=funnode(fspace);

for i=1:m
  tic;
  c=finsolve(model,fspace,stype{i},s,N(i));
  tl(i)=toc;
  C(:,i)=c(:,end);
end
VL=funeval(C,fspace,S);

ts=zeros(1,m);
N=[1 250 75 75 75];

fspace=fundefn('spli',n,0,2*K);    
s=funnode(fspace);
for i=1:m
  tic;
  c=finsolve(model,fspace,stype{i},s,N(i));
  ts(i)=toc;
  C(:,i)=c(:,end);
end
VS=funeval(C,fspace,S);

disp('Maximum errors on [0,2K] for lines, explicit, implicit, CN and stiff methods')
disp('Piecewise Linear')
show(max(abs(premium(:,ones(5,1))-VL)),6)
disp('Cubic Spline')
show(max(abs(premium(:,ones(5,1))-VS)),6)

disp('Timing for lines, explicit, implicit, CN and stiff methods')
disp('Piecewise Linear')
show(tl)
disp('Cubic Spline')
show(ts)

prtfigs(mfilename,'Black-Scholes Option Pricing Model: Approximation Errors',[2])

%% 11.3

% Define parameters
r     = 0.05;     % risk free interest rate
delta = 0;        % rate of dividend payments
kappa = 1;        % mean-reversion parameter on volatility 
m     = 0.05;     % long-run mean volatility
sigma = 0.2;      % volatility of volatility
rho   = -0.5;     % correlation between price and volatility
K     = 1;        % strike price
put   = 1;        % 0=call, 1=put
T     = 1;        % time to maturity

% Create model variable
clear model
model.func='mffin03';
model.T=T;
model.american=0;
model.params={r,delta,kappa,m,sigma,rho,K,put};

% Define approximation space
n=[100 10];
smin=[log(0.01*K) 0.1*m];
smax=[log(5*K) 4*m];
fspace=fundefn('spli',n,smin,smax);
s=funnode(fspace);

% Call solution algorithm
N=fix(T*365+1);            % use one day time steps
N=1000;
c=finsolve(model,fspace,'implicit',s,N);

% Compute BS solution and produce plots
p=linspace(log(.5*K),log(1.5*K),251)';
nu=m;
S=exp(p);

V1=funeval(c(:,end),fspace,{p,nu});
V2=bs(sqrt(nu),S,K,r,delta,T,put);

close all
figure(3);
plot(S,V1,S,V2)
title('Option Values')
xlabel('S')
ylabel('Premium')
xlim([.75*K,1.5*K]);
legend('SV','Black-Scholes')

isigma=impvol(V1,S,K,r,delta,T,put,0);
figure(3);
plot(1./S,isigma,1./S,sqrt(nu)+zeros(size(S,1),1))
title('Implied Volatilities')
xlabel('K')
ylabel('Volatility')
xlim([0.75*K,1.5*K])
legend('SV','Black-Scholes')

nu=linspace(.1*m,4*m,7)';
nu=[.005 .05 .1 .125 .15 .175 .2]';
S=gridmake(p,nu);
V=reshape(funeval(c(:,end),fspace,S),251,7);
figure(3)
plot(exp(p),V)
title('Option Values for Alternative Values of \nu')
xlabel('S')
ylabel('Premium')
nn=length(nu);
legend([repmat('\nu = ',nn,1) reshape(sprintf('%5.3f',nu),5,nn)'],1)

prtfigs(mfilename,'Solution of the Stochastic Volatility Option Pricing Model',[1 2 3])

%% 11.4

% Define parameters
r      = 0.05;               % risk free interest rate
deltaS = 0;                  % dividend rate 
sigma  = 0.2;                % volatility  
K      = 1;                  % exercise price 
put    = 1;                  % put indicator
T      = 1;                  % time-to-maturity

% Create model variable
clear model
model.func='mffin04';
model.T=T;
model.american=1;
model.params={r,deltaS,sigma,K,put};

% Define approximation space
n=100; 
fspace=fundefn('lin',n,0,2*K);    
s=funnode(fspace);

% Call solution algorithm
N=75;
optset('finsolve','keepall',1);
c=finsolve(model,fspace,[],s,N);

% Compute Barone-Adesi/Whaley solution and plot solution and errors
S=linspace(fspace.a,fspace.b,301)';
V=funeval(c(:,end),fspace,S);
[Vbaw,bs,sstarB]=baw(sigma,S,K,r,deltaS,T,put);

figure(4)
plot(S,V,'k-',S,bs,'r-',sstarB,funeval(c(:,end),fspace,sstarB),'k*');
title('Option Premium')
xlabel('S')
ylabel('Premium')
legend('American','European')

figure(4)
plot(S,[Vbaw-V],'k')
title('Approximation Error')
xlabel('S')
ylabel('Error')

% Compute and plot the optimal exercise boundary
V0=feval(model.func,'V0',s,model.params{:});
temp=funeval(c,fspace,s)==V0(:,ones(1,N+1)) & s(:,ones(1,N+1))<=K;
% use upper bound on S* at time values where the upper bound changes
sstar=s(sum(temp)'+1);
sstarl=s(sum(temp)');
tau=linspace(0,T,N+1)';
ind=[1;find(diff(sstar)~=0)+1;N+1];
sstari=sstar(ind);
taui=tau(ind);
% end point adjustments
sstari(1)=K;
send=sstari(end-1)+(sstari(end-1)-sstari(end-2))/(taui(end-1)-taui(end-2))*(taui(end)-taui(end-1));
sstari(end)=(sstari(end)+send)/2;

figure(4)
plot(taui,sstari,'k')
hold on; stairs(tau,sstar,'r--'); stairs(tau,sstarl,'r--'); hold off
title('Early Exercise Boundary')
xlabel('\tau')
ylabel('S^*')
ylim([0.8 1.05])

prtfigs(mfilename,'Solution of the American Put Option Pricing Model',[1 2 3])

%% 11.5

% Define parameters
r      = 0.05;               % risk free interest rate
delta  = 0;                  % dividend rate 
sigma  = 0.2;                % volatility  
K      = 1;                  % exercise price 
put    = 0;                  % put indicator
T      = 1;                  % time-to-maturity
Sb     = 0.8;                % barrier 

% Create model variable
clear model
model.func='mffin04';
model.T=T;
model.american=0;
model.params={r,delta,sigma,K,put};
model.barrier=[0 Sb 1;0 inf 1];

% Define approximation space
n=75; 
fspace=fundefn('spli',n,Sb,2*K);    
S=funnode(fspace);

% Call solution algorithm
[c,Vhat]=finsolve(model,fspace);

% Create plots
S=sort([linspace(.5*K,2*K,501)';Sb-eps;Sb]);
Vhat=funeval(c,fspace,S);
Vhat(S<Sb)=0;
Vbs=bs(sigma,S,K,r,delta,T,put);
V=Vbs-(S./Sb).^(1-2*r/sigma^2).*bs(sigma,Sb^2./S,K,r,delta,T,put);
V(S<Sb)=0;

figure(5)
plot(S,Vhat,S,Vbs);
title('Down-and-out Call Option Premium')
xlabel('S')
ylabel('Pemium')
legend({'Barrier Option','Vanilla Option'},2);
xlim([.5 1.5]*K)
ylim([0 0.5])

figure(5)
plot(S,[V-Vhat])
title('Approximation Error')
xlabel('S')
ylabel('Error')
xlim([Sb,2*K])

disp('Sb Vb(Sb) Vbs(Sb)')
disp([S(S==Sb) Vhat(S==Sb) Vbs(S==Sb)])

n=75; 
lambda=.5;
Delta=(2*K-Sb)/(n-2);
fspace=fundefn('lin',n,Sb-lambda*Delta,2*K);    
% Call solution algorithm
s=funnode(fspace);
s(abs(s-Sb)<1e-13)=Sb;
[c,Vhat]=finsolve(model,fspace,[],s);

Vhat=funeval(c,fspace,S);
figure(5)
plot(S,[V-Vhat])
title('Approximation Error with Barrier Not a Node')
xlabel('S')
ylabel('Error')
xlim([Sb,2*K])

prtfigs(mfilename,'Solution of the Barrier Option Pricing Model',[1 2 3])

%% 11.6

% Define parameters 
kappa = 0.1;       % speed of mean reversion
alpha = 0.05;      % long-run mean interest rate
sigma = 0.1;       % interest rate volatility 
TB    = 30;        % bond maturity
K     = 0.2;       % exercise price
put   = 0;         % put indicator
TO    = 1;         % option maturity

% Create model variable for the bond
clear model
modelB.func='mffin06';
modelB.T=TB;
modelB.american=0;
modelB.params={kappa,alpha,sigma};

% Define approximation space for the bond
n=20;
fspaceB=fundefn('cheb',n,0,2);

% Call the solver
cB=finsolve(modelB,fspaceB);

% Create model variable the option
modelO.func='mffin06';
modelO.T=TO;
modelO.american=0;
modelO.params={kappa,alpha,sigma,K,put,cB,fspaceB};

% Define approximation space for the option
n=80;
fspaceO=fundefn('spli',n,0,2);

% Call the solver
cO=finsolve(modelO,fspaceO);

% Create plots
x=linspace(0,.25,201)';
Vhat=funeval(cO,fspaceO,x);

figure(6)
plot(x,Vhat); 
title(['As a Function of Short Rate'])
xlabel('r')
ylabel('Option Premium')
ylim([0 0.25])

BondVal=funeval(cB,fspaceB,x);
figure(6)
plot(BondVal,Vhat); 
title(['As a Function of Bond Value'])
xlabel('Bond Value')
ylabel('Option Premium')
ylim([0 0.25])
xlim([0 .45])

prtfigs(mfilename,'Solution of the Compound Option Pricing Model',[1 2])

%% 11.7

% Define parameters
r     = 0.1;
delta = 0;
sigma = 0.15;
L     = 1/2;
put   = 1;
tau   = 1/4;


tau=1/4;


% Create model variable
clear model
model.func='mffin07';
model.T=tau;
model.american=0;
model.params={r,delta,sigma,L,put};

% Define approximation space
n=201; 
fspace=fundefn('lin',n,0,1);

% Call solution algorithm
c=finsolve(model,fspace);

M=[.5 .75 1 1.25 1.5];
funeval(c,fspace,(L-tau)*M')'
return
% Transform solution to natural state space (y to S)
S=linspace(0.001,2,101)';
M=[.5 .75 1 1.25 1.5];
m=length(M);
Vhat=zeros(length(S),m);
for i=1:m
  y=((L-tau)*M(i))./S;
  Vhat(:,i)=S.*funeval(c,fspace,y);
  if ~put, Vhat(y>fspace.b,i)=0; end
end

% Create plots
figure(7)
plot(S,Vhat);
title('Call Option Premium')
xlabel('S')
ylabel(['V(S,M,' num2str(tau) ')']);
nn=length(M);
legend([repmat('M = ',nn,1) reshape(sprintf('%4.2f',M),4,nn)'],2)
set(gca,'ylim',max(0,get(gca,'ylim')));

y=linspace(fspace.a,fspace.b,101)';
figure(7)
plot(y,funeval(c,fspace,y));
title('Approximation Function')
xlabel('y')
ylabel('v(y,\tau)');
ylim([-0.1 0.6])

%% 11.8

% Define parameters
T     = 30;         % time-to-maturity
kappa = 0.1;        % speed of mean reversion
alpha = 0.05;       % long-run mean interest rate 
sigma = 0.1;        % interest rate volatility

% Convert parameters to standard form
a=kappa*alpha;
A=-kappa;
C=sigma;
b=0;
B=1;

% Call affine asset pricing solver
tau=linspace(0,T,301)';
[beta,beta0]=affasset(tau,a,A,b,B,C,1,0,0,0);

% Create plots
r=linspace(0,0.25,101)';
V=exp(beta0(end)+beta(end)*r);

figure(8)
plot(r,V); 
title([num2str(T) ' Year Zero-Coupon Bond Price'])
xlabel('r')
ylabel('V(r)')

figure(8)
plot(tau,[beta beta0]); 
title('\beta and \beta_0')
xlabel('Time to Maturity')

r=(0.03:0.01:0.08);
m=length(r);

warning off
R=-(beta0(:,ones(m,1))+beta*r)./tau(:,ones(m,1));
warning on
R(1,:)=r;
figure(3)
plot(tau,R)
xlabel('Time to Maturity')
ylabel('Yield')
title('Term Structures for Alternative Short Rates')

nn=length(r);
legend([repmat('r = ',nn,1) reshape(sprintf('%4.2f',r),4,nn)'])

prtfigs(mfilename,'Term Structures for Alternative Short Rates',[3])

%% 11.9

% Define parameters
kappa = 0.0363;
alpha = 0.0692; 
sigma = 0.0272;
 
% Create model variable
clear model
model.func='mffin01';
model.T=30;
model.american=0;
model.params={kappa,alpha,sigma};

% Define approximation space
n=20;
fspace=fundefn('cheb',n,0,2);
s=funnode(fspace);

% Call solution algorithm
optset('finsolve','keepall',1);
c=finsolve(model,fspace,'lines',s,120);

% Calibrate to the data
y=[4.44	 4.49 	4.51	4.63	4.63	4.62	4.82	4.77	5.23];
tau=[.25 .5 1 2 3 5 7 10 30];
V=exp(-y/100.*tau);
m=length(V);
t=(0:0.25:30);
tind=[2 3 5 9 13 21 29 41 121];
% t(tind)=tau and columns of c correspond to t
s=findstate(c(:,tind),fspace,alpha,V);

% Create plots
Vhat=funeval(c,fspace,s);
warning off
yhat=-100*log(Vhat)./t; yhat(1)=100*s;
warning on

figure(9)
plot(tau,y,'*',t,yhat)
title('Actual and Fitted Bond Yields')
xlabel('Time to Maturity')
ylabel('Yield')

disp('Model short interest rate and 3-month rate')
disp([s*100 y(1)])

prtfigs(mfilename,'Actual and Fitted Bond Yields',[1])

%% 11.10

% Define problem parameters
  alpha = 0.14;
  delta = 0.02;
  gamma = 0.5;
  rho   = 0.05;

% Pack model structure
  model.func='mfsc01';
  model.params={alpha,delta,gamma,rho};
  
% Define nodes and basis matrices
  n=20;
  smin=0.2; 
  smax=2;
  fspace=fundefn('cheb',n,smin,smax);
  snodes=funnode(fspace);

% Define initial values
  v0=((rho*snodes).^(1-gamma)-1)/(1-gamma)/rho;
  x0=(rho+delta)*snodes;
  
% Call solver  
  [cv,s,v,x,resid]=scsolve(model,fspace,snodes,v0);
  
% Get steady state
  [Kstar,Cstar]=ctsteadystate(model,fspace,cv);

% Display steady state results
  disp('Steady State Capital and Consumption')
  disp([Kstar Cstar])

% Plot value function
  figure(10)
  plot(s,v)
  title('Value Function')
  xlabel('K')
  ylabel('V(K)')
  xlim([0 2])

% Plot optimal control function
  figure(10control)
  plot(s,x,Kstar,Cstar,'k*')
  title('Optimal Consumption Rule')
  xlabel('K')
  ylabel('C')
  xlim([0 2])

% Plot residual function
  figure(10resid)
  plot(s,resid)
  title('Approximation Residual')
  xlabel('K')
  ylabel('Residual')
  xlim([0 2])
  