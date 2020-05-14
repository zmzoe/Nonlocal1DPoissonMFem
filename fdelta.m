
format long

k=1;
delta=0.000001;

ureal = @(x) cos(2*pi*k*x); 
f=@(x) -(2*cos(2*pi*k*x)*(cos(2*pi*delta*k) - 2*delta^2*k^2*pi^2 + 2*delta*k*pi*sin(2*pi*delta*k) - 1))/(delta^4*k^2*pi^2);

g=@(x) 4*pi^2*cos(2*pi*x);

 
h=@(x) f(x)-g(x);

