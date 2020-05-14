
syms x s delta

G1=  (-2/delta^2)*x^2 + (4/delta)*x - (2*(delta/4 + (delta*(2*delta - 1))/4))/delta^2;
G2=  1;
G3=  (2/delta^2)*x^2 + (-(2*(2*delta + 1))/delta^2)*x + (2*((3*delta)/4 + (delta*(2*delta + 1))/4 + 1/4))/delta^2;
G4=  -1;


%% if (0,delta)
J12=G1*delta;

%if delta>1/4
J11a=int(subs(G1,x,x+s),s,0,delta-x)+int(subs(G2,x,x+s),s,delta-x,1/2-x)+int(subs(G3,x,x+s),s,1/2-x,delta);

F1a=-2/(delta^2)*(J11a-J12);

%if delta<1/4
J11b=int(subs(G1,x,x+s),s,0,delta-x)+int(subs(G2,x,x+s),s,delta-x,delta);

F1b=-2/(delta^2)*(J11b-J12);

%% if (delta,1/2)
J22=G2*delta;

J21=int(subs(G2,x,x+s),s,0,1/2-x)+int(subs(G3,x,x+s),s,1/2-x,delta);

F2=-2/(delta^2)*(J21-J22);

%% if (1/2, delta+1/2)
J32=G3*delta;

%if delta>1/4
G11=subs(G1,x,x-1);
J31a=int(subs(G3,x,x+s),s,0,1/2+delta-x)+int(subs(G4,x,x+s),s,1/2+delta-x,1-x)+int(subs(G11,x,x+s),s,1-x,delta);

F3a=-2/(delta^2)*(J31a-J32);

%if delta<1/4
J31b=int(subs(G3,x,x+s),s,0,1/2+delta-x)+int(subs(G4,x,x+s),s,1/2+delta-x,delta);

F3b=-2/(delta^2)*(J31b-J32);


%% if (delta+1/2,1)
J42=G4*delta;
J41=int(subs(G4,x,x+s),s,0,1-x)+int(subs(G11,x,x+s),s,1-x,delta);

F4=-2/(delta^2)*(J41-J42);




