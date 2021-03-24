clc;
clear;
t = zeros(1,100);
for n = 1:200
    A = rand(n,n);
    b = rand(n,1);
    tic;
    x = A\b;
    t(1,n) = toc;
    N(n)=n;
end
[N' (t*1000)']
plot(t)


n=4000
A = rand(n,n);
b = rand(n,1);
t1=tic;
x = A\b;
t2=toc
   