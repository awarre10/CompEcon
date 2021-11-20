function [v,x,vjac] = mymax(s,x,c,e,maxit,tol,fspace,w,alpha,beta,gamma,delta)
[xl,xu] = myfunc('b',s,x,e,alpha,beta,gamma); 5bfunc(s);
K = length(e);
for it=1:maxit
    [f,fx,fxx] myfunc('f',s,x,e,alpha,beta,gamma);
    Ev = 0; Evx = 0; Evxx = 0;
    for k = 1:K
        [g,gx,gxx] = myfunc('g',s,x,e(k),alpha,beta,gamma);
        vn = funeval(c,fspace,g);
        vnder1 = funeval(c,fspace,g,1);
        vnder2 = funeval(c,fspace,g,2);
        Ev = Ev + w(k)*vn;
        Evx = Evx + w(k)*(vnder1.*gxx + vnder2.*gx.^2);
    end
    v = f+delta*Ev;
    delx =-(fx+delta*Evx)./(fxx+delta*Evxx);
    delx = min(max(delx,xl-x),xu-x);
    x = x + delx;
    if norm(delx)<tol, break, end
end
vjac 0
for k=1:K
    g = myfunc('g',x,x,e(k),alpha,beta,gamma);
    vjac = vjac + delta*x(k)*funbas(fspace,g);
end
