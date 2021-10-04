function main
ex_1()
ex_2()
ex_3()
ex_4()
ex_5()
end

% Excercise 1 %

function ex_1
% a:
syms p
qd = 2*p^(-0.5);
ex_1a = int(qd,1,4)

% b:
n= 11; a= 1; b= 4;
[x,w] = qnwtrap(n,a,b);
trap = w'*(qd_i(x));

disp("(b) numerically using an 11-node trapezoid rule.")
disp(trap)
% c:
[x,w] = qnwsimp(n,a,b);
simp = w'*(qd_i(x));

disp("(c) numerically using an 11-node Simpson rule.")
disp(simp)

% d:
[x,w] = qnwlege(n,a,b);
lege = w'*(qd_i(x));

disp("(d) numerically using an 11-node Gauss-Legendre rule.")
disp(lege)

% e:
[x,w] = qnwequi(n,a,b, 'N'); % Neiderriert rule
equi = w'*(qd_i(x));

disp("(e) numerically using an 11-node equidistributed sequence rule.")
disp(equi)

end

function [y] = qd_i(p)
y = 2*p.^(-0.5);
end

% Exercise 2 %

function ex_2
clear all; clc;
disp("Write a program that solves numerically the following expression for α:");
alpha = bisect('fun_ction',0,1)
end

function lambda=fun_ction(alpha)
[lambda, w] = qnwtrap(20000, 0, 10000);
trap = sum(w.*exp(alpha*lambda-lambda.^2./2));
lambda = alpha*trap-1;
end

% Exercise 3 %

function ex_3()
% Monte Calro Method
n = 1000;
x = -log(1-rand(n,1)); % F^(-1) = -log(1-U) F(x)= 1-exp(x)
disp(mean(x.^2))

n = 10000;
x = -log(1-rand(n,1));
disp(mean(x.^2))

n = 100000;
x = -log(1-rand(n,1));
disp(mean(x.^2))
%{
The expectation of a function of x is roughly equal to (1/n)sum(from 1 to
n)f(xi)^2. To use this method go through the inverse transform.
Inverse CDF: f(x)
f(x)^(-1)
%}

% Neiderreiter quasi-Monte Carlo integration
n = [1000; 10000; 100000];
% f(x) = 1-exp(-x)
% f(x) = (-log(1-x))
disp("Neiderreiter quasi-Monte Carlo integration")
for i = 1:length(n)
    [x,~] = qnwequi(n(i),0,1,'N');
    mean((-log(1-x)).^2)
end
end

function [y3] = fun3(x)
y3 = 1-exp(-x);
end

% Exercise 4 %

function ex_4()
%   The centered finite difference method is more accurate because its
%   error is one order more accurate than the finite one sided difference.
%{
a = 0; b = 1/2h; c = -1/2h;
a*(f) + (f(x+h))/2h - (f(x+lambda*h))/2h = fprime(x) + o*h^2
%}
fj = fdjac('f4', 1);
fd = fdhess('f4',1);
end

% Exercise 5 %

function ex_5()
close all

disp(' ')
disp('DEMDIF03 Commercial Fishery Model (from V.L. Smith)')
alpha = 1; % Parameter values
% beta=2.75; f=0.025; del=10;            

% Trajectories
maxt=30;                                 
n=100;   
t=[0:0.025:4.95 5:(maxt-5)/n:maxt]';
xvals=0:0.2:2;
yvals=0:0.2:2;
x0=phase1(xvals,yvals);                  
x0=x0(:,find(x0(1,:)>0));               % eliminate values with S=0
    % time values
[t,x]=rk4('lv',t,x0,[],alpha); % call the ODE solver

% Isoclines
one = 1;
s=0.01:0.01:1;
x1=(one)./(alpha-s);            % x one-isocline
x2=(one)./(s-one);         % y one-isocline

% Plot the phase diagram
figure(1)
plot(squeeze(x(:,1,:)),squeeze(x(:,2,:)),'-')
hold on
plot(s,x1,':');
plot(s,x2,':');
hold off
xlabel('x');
ylabel('y');
title('Phase Diagram for Predator-Prey Model');
axis([0 2 0 2]);

h=text(.005,.3,'K''<0');
set(h,'fontsize',7)
h=text(.06,.1,'K''>0');
set(h,'fontsize',7)
h=text(.825,.15,'S''>0');
set(h,'fontsize',7)
h=text(.825,.35,'S''<0');
set(h,'fontsize',7)

h=text(0.085,1.125,'A');
set(h,'fontsize',8)
h=text(.325,1.8,'B');
set(h,'fontsize',8)
h=text(0.5,1.5,'C');
set(h,'fontsize',8)

prtfigs(mfilename,'Phase Diagram for Predator-Prey Model',1)



% PHASE   A utility to generate boundary points for 2D phase plot.
% INPUTS: X,Y vectors of values for the x and y axes.
% OUTPUT: P, a 2xK matrix of initial values that can be passed to
%   an ODE solver (e.g., RK4).


end

function dk = lv(t,k,flag,alpha)
x=k(1,:);
y=k(2,:);
dx=alpha*x-x.*y;
dy=x.*y-y;
dk=[dx;dy];
end

function p=phase1(x,y);
  miny=min(y);
  maxy=max(y);
  x=x(:)';
  y=y(:)';
  nx=length(x);
  ny=length(y);
  y=y(2:ny-1);
  ny=ny-2;
  p=[x;(miny+zeros(1,nx))];
  if ny>0 p=[p [(max(x)+zeros(1,ny));y]]; end
  p=[p [fliplr(x);(maxy+zeros(1,nx))]];
  if ny>0 p=[p [(min(x)+zeros(1,ny));fliplr(y)]]; end
end
% The center is at (1,1),
% Like the top of a mountain on an elevation diagram.
% Just enough prey to sustain the relationship.
function x=bisect(f,a,b,varargin)
% get options
  tol         = optget('bisect','tol',1e-4);
  checks      = optget('bisect','checks',0);
  mustbracket = optget('bisect','mustbracket',1);

% Perform checks
  if checks
    if nargin<3
      error('At least three parameters must be passed to BISECT');
    end
    if size(a)~=size(b)
      error('In BISECT: Lower and upper ranges must be the same size');
    end
    if any(a>b)
      error('Lower bound greater than upper bound');
    end
  end
  sa=sign(feval(f,a,varargin{:}));
  sb=sign(feval(f,b,varargin{:}));
  if any(sa==sb) & mustbracket
    error('In BISECT: root not bracketed')
  end

% Initializations  
  dx = 0.5*(b - a);
  tol = dx.*tol;
  x = a + dx;
  dx = sb.*dx;

% Iteration loop
  while any(abs(dx)>tol)
    dx = 0.5*dx;
    x = x - sign(feval(f,x,varargin{:})).*dx;
  end

    return
end

function H = fdhess(f,x,varargin)
  tol      = optget(mfilename,'tol',eps.^(1/4));
  diagonly = optget(mfilename,'diagonly',0);
 
  k = size(x,1);
  fx = feval(f,x,varargin{:});
 
  % Compute the stepsize (h)
  h = tol*max(abs(x),1);
  xh = x+h;
  h = xh-x;    
  ee = sparse(1:k,1:k,h,k,k);
 
  % Compute forward and backward steps
  gplus = zeros(k,1);
  gminus = zeros(k,1);
  for i=1:k
    gplus(i) = feval(f,x+ee(:,i),varargin{:});
    gminus(i) = feval(f,x-ee(:,i),varargin{:});
  end
   
  H=h*h';
  % Compute double steps
  if diagonly
    for i=1:k
      H(i,i) = (gplus(i)+gminus(i)-2*fx)/ H(i,i);
    end
  else
    for i=1:k
      for j=1:k
        if i==j
          H(i,j) = (gplus(i)+gminus(j)-2*fx)/ H(i,j);
        else
          fxx=feval(f,x+ee(:,i)-ee(:,j),varargin{:});
          H(i,j) = (gplus(i)+gminus(j)-fx-fxx)/ H(i,j);
        end
      end
    end
    H=(H+H')/2;
  end
end
function fjac = fdjac(f,x,varargin)

tol    = optget(mfilename,'tol',eps.^(1/3));

h = tol.*max(abs(x),1);
xh1=x+h; xh0=x-h;
h=xh1-xh0;
for j=1:length(x);
   xx = x;
   xx(j) = xh1(j); f1=feval(f,xx,varargin{:});
   xx(j) = xh0(j); f0=feval(f,xx,varargin{:});
   fjac(:,j) = (f1-f0)/h(j);
end
end
