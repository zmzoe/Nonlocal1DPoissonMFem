%h=delta/M
%Omega=(0,1) periodical boundary condition

%delta=0.5
%M=1
%h=0.5
%#points=3



syms x x_1 x_2 x_3 x_4 p_1 p_2 d h s v_1 v_2



phi_1minus=(x_2-x)/h;

phi_2plus=(x-x_1)/h;
phi_2minus=(x_3-x)/h;

phi_3plus=(x-x_2)/h;
phi_3minus=(x_4-x)/h;

phi_4plus=(x-x_3)/h;


%% ~~~~~~~~~~x1 < x < x2~~~~~~~~~~~~~~~G_delta^* . d=delta

g11=p_1*(subs(phi_1minus,x,x+s))+p_2*(subs(phi_2plus,x,x+s)); % 0 < s < x_2-x
g12=p_2*(subs(phi_2minus,x,x+s))+p_1*(subs(phi_3plus,x,x+s)); % x_2-x < s < delta

f11=2/(d^2)*int(g11,s,0,x_2-x);
f12=2/(d^2)*int(g12,s,x_2-x,d);

g13=p_1*phi_1minus+p_2*phi_2plus;

f13=2/(d^2)*int(g13,s,0,d);%remark minus f13 

f1=f11+f12-f13; 

%collect(f1,p_1)


%% ~~~~~~~~~~x2 < x < x3~~~~~~~~~~~~~~~G_delta^* . d=delta

g21=p_2*(subs(phi_2minus,x,x+s))+p_1*(subs(phi_3plus,x,x+s)); % 0 < s < x_3-x
g22=p_1*(subs(phi_3minus,x,x+s))+p_2*(subs(phi_4plus,x,x+s)); % x_3-x < s < delta

f21=2/(d^2)*int(g21,s,0,x_3-x);
f22=2/(d^2)*int(g22,s,x_3-x,d);

g23=p_2*phi_2minus+p_1*phi_3plus;

f23=2/(d^2)*int(g23,s,0,d);

f2=f21+f22-f23; 

%collect(f2,p_1)

%% ~~~~~~~~~~~~~~~Integration~~~~~~~~~~~~~~~~~~~~~~~
%\|Gp\|^2
I11=int(f1^2,x,x_1,x_2);
I12=int(f2^2,x,x_2,x_3);

L11=I11+I12;
L11=subs(L11,[x_1,x_2,x_3,x_4],[0,h,2*h,3*h]);

%\|p\|^2
I21=int(g13^2,x,x_1,x_2);
I22=int(g23^2,x,x_2,x_3);

L12=I21+I22;
L12=subs(L12,[x_1,x_2,x_3,x_4],[0,h,2*h,3*h]);

L1=L11+L12;
%collect(L1,p_1)

%% frac_up B
I21=int(f1,x,x_1,x_2);
I22=int(f2,x,x_2,x_3);

L2=v_1*I21+v_2*I22;
L2=subs(L2,[x_1,x_2,x_3,x_4],[0,h,2*h,3*h]);
collect(L2,v_1)

