% PS#2

% Exercise 3
stp=0.02; strt=-2; fnsh=3;
x = strt:stp:fnsh;
f = 1-exp(3*x)
A = [0 -1 2;-2 -1 4;2 7 -3]
B = [-7 1 1;7 -3 -2;3 5 0]
y = [3;-1;2]
C = A*B
C\y
C = A.*B
C\y

% Exercise 4
t = 1970:2025
e = 0.2*randn(size(t))
y = 5+0.05*t+e
plot(t,y)

% Exercise 5
% a:
E = 0.5*0.7+0.5*1.3
V = 0.5*(0.7-1)^2+0.5*(1.3-1)^2
% b:
y = randsrc(2,1,[0.7 1.3; 0.5 0.5])
w = 0.5*ones(size(y))
a = 1
for it=1:100;
    aold=a;
    p = 3-2*a*y;
    a = 0.5+0.5*f;
    if abs(a-aold)<1.e-8, break, end
end
disp(a); disp(w'*p); disp(f)
V=0.5*[max(0.7-1)-f]^2+0.5*[max(1.3-1)-f]^2
