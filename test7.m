syms x delta k

ureal = @(x) sin(2*pi*x);
du=diff(ureal,x,1);
f=-diff(du,x,1);