
% Exercise 1

% a:

% f'(x) = 2x-1

% b: 

function y = j(x)
    y = (x^2)-c;
    x = bisect('j',1,2);
end

% c:

function x = newtroot(c)
x-1;
tol = 5E-5;
maxit = 50;
for it = 1:50
    fval = (x^2)-c
    fjac = 2x
    x=x-fjac\fval;
    if norm(fval)<tol,break,end
end
    
% Exercise 2

% a:

function V = BSVal(S,K,tau,r,delta,sigma)
d = (log(exp(-delta*tau)*S)-log(exp(-r*tau)*K))/(sigma*sqrt(tau)) + (1/2)*sigma*sqrt*(tau);
V = exp(-delta*tau)*S*cdfn(d) - exp(-r*tau)*K*cdfn(d-sigma*sqrt(tau));

% b:

function sigma = ImpVol(S,K,tau,r,delta,V)
sigma = 0.1 % initial guess for volatility
S = 60;
K = 58; 
tau = 0.5; 
r = 0.045;
delta = 0.0125 % divident yield percentage
V = 4.77;
tol = 5E-5;
maximt = 50;
for it=1:maxit
    fval = V - BSVal(S,K,tau,r,delta,sigma)
    d = (log(exp(-delta*tau)*sqrt(tau/(2*pi))*exp(-0.5*d^2);
    sigma = sigma - fjac\fval;
    if coder.extrinsic norm(fval) < tol, break, end
end


% Exercise 3

function [fval,fjac] = f3(x)
    fval = [200x(1)(x(2)-x(1)^2)-x(1)+1,100((x(1)^2)-x(2))];
    fjac = [[-600*(x(1)^2)+200*x(2)-1, 200*x(2); 200*x(1), -100];
end
end
newton('f3'[1.25,13])
broyden('f3'[1.25,13])
end
