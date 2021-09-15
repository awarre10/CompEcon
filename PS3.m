% Excercise 1

% a:

A = [54 14 -11 2; 14 50 -4 29; -11 -4 55 22; 2 29 22 95];
B = [1; 1; 1; 1];
[L,U] = lu(A);
y = L\B
x = U\y

% b:

x_GJ = gjacobi(A,B);

% c:

x_GS = gseidel(A,B);

% Exercise 2

% a:

A = randn(100,100);
b = randn(100,1);

tic
x=A\b;
toc
tic
for i=1:10
    x=A\b;
end
toc
tic
for i=1:50
    x=A\b;
end
toc

% b:

tic
[L,U] = lu(A);
x = U\(L\b)
toc
tic
for i=1:10
    x = U\(L\b);
end
toc
tic
for i=1:50
    x = U\(L\b);
end
toc

% c:

tic
c=(inv(A))*b;
toc
tic
for i=1:10
    c;
end
toc
tic
for i=1:50
    c;
end
toc