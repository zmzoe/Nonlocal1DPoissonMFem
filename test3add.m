syms x s delta 

u1= -x+1/4;
u2= x-3/4;

G1=2/(delta^2)*(int(subs(u1,x,x-s),s,0,x)+int(subs(u2,x,x+1-s),s,x,delta)-u1*delta);
G2=2/(delta^2)*(int(subs(u1,x,x-s),s,0,delta)-u1*delta);
G3=2/(delta^2)*(int(subs(u2,x,x-s),s,0,x-1/2)+int(subs(u1,x,x-s),s,x-1/2,delta)-u2*delta);
G4=2/(delta^2)*(int(subs(u2,x,x-s),s,0,delta)-u2*delta);

%delta>1/4 
f1=-2/(delta^2)*(int(subs(G1,x,x+s),s,0,delta-x)+int(subs(G2,x,x+s),s,delta-x,1/2-x)+int(subs(G3,x,x+s),s,1/2-x,delta)-G1*delta);
f2=-2/(delta^2)*(int(subs(G2,x,x+s),s,0,1/2-x)+int(subs(G3,x,x+s),s,1/2-x,delta)-G2*delta);
f3=-2/(delta^2)*(int(subs(G3,x,x+s),s,0,1/2+delta-x)+int(subs(G4,x,x+s),s,1/2+delta-x,1-x)+int(subs(G1,x,x+s-1),s,1-x,delta)-G3*delta);
f4=-2/(delta^2)*(int(subs(G4,x,x+s),s,0,1-x)+int(subs(G1,x,x+s-1),s,1-x,delta)-G4*delta);

%delta<=1/4 ĞèÒª·ÖÁù¶Î
%(0,delta) 
%(delta,1/2-delta) 
%(1/2-delta,1/2) 
%(1/2,1/2+delta) 
%(1/2+delta,1-delta) 
%(1-delta,1) 

ff1=-2/(delta^2)*(int(subs(G1,x,x+s),s,0,delta-x)+int(subs(G2,x,x+s),s,delta-x,delta)-G1*delta);
ff2=-2/(delta^2)*(int(subs(G2,x,x+s),s,0,delta)-G2*delta);
ff3=-2/(delta^2)*(int(subs(G2,x,x+s),s,0,1/2-x)+int(subs(G3,x,x+s),s,1/2-x,delta)-G2*delta);
ff4=-2/(delta^2)*(int(subs(G3,x,x+s),s,0,1/2+delta-x)+int(subs(G4,x,x+s),s,1/2+delta-x,delta)-G3*delta);
ff5=-2/(delta^2)*(int(subs(G4,x,x+s),s,0,delta)-G4*delta);
ff6=-2/(delta^2)*(int(subs(G4,x,x+s),s,0,1-x)+int(subs(G1,x,x+s-1),s,1-x,delta)-G4*delta);











