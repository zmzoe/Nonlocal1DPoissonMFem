
syms x delta s 

u1= -x+1/4;
u2= x-3/4;

%% if [0,1/2]
I12=(-x+1/4)*delta;
% if [0,delta]
u22=subs(u2,x,x+1);%move to left
I11a=int(subs(u1,x,x-s),s,0,x)+int(subs(u22,x,x-s),s,x,delta);

G1=2/(delta^2)*(I11a-I12);
% if [delta,1/2]
I11b=int(subs(u1,x,x-s),s,0,delta);
G2=2/(delta^2)*(I11b-I12);

%% if [1/2,1]
I22=(x-3/4)*delta;
% if [1/2, 1/2+delta]
u11=subs(u1,x,x-1);%move to left
I21a=int(subs(u2,x,x-s),s,0,x-1/2)+int(subs(u11,x,x-s),s,x-1/2,delta);

G3=2/(delta^2)*(I21a-I22);
% if [delta+1/2,1]
I21b=int(subs(u2,x,x-s),s,0,delta);
G4=2/(delta^2)*(I21b-I22);



