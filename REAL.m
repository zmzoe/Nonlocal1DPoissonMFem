
%different real solution 
% G :original definition

%%
syms x s delta k

% k is an interger
u=cos(2*pi*k*x);
p=-2*k*pi*sin(2*pi*k*x);

%% f

Gu=2/(delta^2)*int(subs(u,x,x-s)-u,s,0,delta);

f=-2/(delta^2)*int(subs(Gu,x,x+s)-Gu,s,0,delta);

%% g

g=p-Gu;

%% HH
HH=-2/(delta^2)*int(subs(p,x,x+s)-p,s,0,delta);


