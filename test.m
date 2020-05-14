



syms x k delta s

%{
f1=@(x)-(x-1/4).^2+1/16;
f2=@(x) (x+1/4).^2-1/16;

u=@(x) subs(f1(x),x,x+1).*(x>=-1 & x<-1/2)+f2(x).*(x>=-1/2 & x<0)+f1(x).*(x>=0 & x<1/2)+subs(f2(x),x,x-1).*(x>=1/2 & x<1)+subs(f1(x),x,x-1).*(x>=1 & x<3/2)+subs(f2(x),x,x-2).*(x>=3/2 & x<=2);
%u=@(x) cos(2*pi*k*x);
plocal=@(x) diff(u(x),x,1);

Gu=@(x) 2/(delta^2)*int(subs(u(x),x,x-s)-u(x),s,0,delta);

f=@(x) -2/(delta^2)*int(subs(Gu(x),x,x+s)-Gu(x),s,0,delta);

g=@(x) plocal(x)-Gu(x);

%HH=@(x) f(x)-2/(delta^2)*int(subs(g(x),x,x+s)-g(x),s,0,delta);

HH=@(x) -2/(delta^2)*int(subs(plocal(x),x,x+s)-plocal(x),s,0,delta);
%}

u=@(x) sin(6*pi*x);   
Gu=@(x) 2/(delta^2)*int(subs(u(x),x,x-s)-u(x),s,0,delta);
f=@(x) -2/(delta^2)*int(subs(Gu(x),x,x+s)-Gu(x),s,0,delta);





