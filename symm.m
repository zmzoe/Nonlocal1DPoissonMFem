%syms

%% Notations

syms d p_i p_i1 p_i2 p_i3 x x_i3 x_i2 x_i1 x_i x_i_1 s h i v_i

phi_iplus=(x-x_i_1)/h;
phi_iminus=(x_i1-x)/h;

phi_ip1plus=(x-x_i)/h;
phi_ip1minus=(x_i2-x)/h;

phi_ip2plus=(x-x_i1)/h;
phi_ip2minus=(x_i3-x)/h;

%%
%                            G_delta^*p(x) and x\in (x_i, x_i+1)
%
%% when x_i<x<x_i1-d  G_d^*p

g1=p_i*(subs(phi_iminus,x,x+s))+p_i1*(subs(phi_ip1plus,x,x+s))-p_i*phi_iminus-p_i1*phi_ip1plus;

f1=2/(d^2)*int(g1,s,0,d);

%collect(f1)

%% when x_i1-d<x<x_i1 G_d^*p

%when 0<s<x_i1-x
g21=g1;
f21=2/(d^2)*int(g21,s,0,x_i1-x);
%collect(f21,(x_i1-x))

%when x_i1-x<s<d
g22=p_i1*(subs(phi_ip1minus,x,x+s))+p_i2*(subs(phi_ip2plus,x,x+s))-p_i*phi_iminus-p_i1*phi_ip1plus;
f22=2/(d^2)*int(g22,s,x_i1-x,d);

f2=f21+f22;

%%
%                            \|G_delta^*p\|_0^2
%
%% Integrate over (x_i, x_i1)

I1=int(f1^2,x,x_i,x_i1-d);
I2=int(f2^2,x,x_i1-d,x_i1);

L_2squar=I1+I2;
L2=subs(L_2squar,[x_i,x_i1,x_i2],[1i*h,(1i+1)*h, (1i+2)*h]);
%collect(L2, p_i);
%%  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~test~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%{
syms delta

a=[(15*h-7*delta)/(15*h^2), (3*delta-5*h)/(5*h^2), -(2*delta)/(15*h^2)];
b=[(3*delta-5*h)/(5*h^2), (15*h-8*delta)/(15*h^2), -delta/(15*h^2)];
c=[-(2*delta)/(15*h^2), -delta/(15*h^2),   delta/(5*h^2)];

H=[a;b;c];
S=[p_i,p_i1,p_i2]*H*[p_i;p_i1;p_i2];
L1=subs(S,delta, d);
collect(L1, p_i)
%}


%%
%                            \|p\|_0^2
%
%% when x_i < x < x_i1    compute p(x)

g3=p_i*phi_iminus+p_i1*phi_ip1plus;


%% Integrate over (x_i, x_i1)
I3=int((g3)^2,x,x_i,x_i1);

L2_P=subs(I3,[x_i,x_i1],[1i*h,(1i+1)*h]);
%collect(L2_P)
%% ~~~~~~~~~~~~~~~~~~test~~~~~~~~~~~~~~~~~~~~~~~~~~~
%{
syms a b p1 p2 c d
H=[p1, p2]*[a b; c d]*[p1;p2];
collect(H)
%}

%%
%                            \|v\|_0^2 v is piecewise constant . EASY
%
%%
%                                   int(G_delta^*p\cdot v)
%
%%
I41=int(f1,x,x_i,x_i1-d); %int_(x_i,x_i1-d)G_d^*p
I42=int(f2,x,x_i1-d,x_i1);%int_(x_i1-d,x_i1)G_d^*p

I4=(I41+I42)*v_i; %v_i piecewise constant
I4=subs(I4,[x_i,x_i1,x_i2],[1i*h,(1i+1)*h, (1i+2)*h]);

L_frac_up=collect(I4);
%collect(L_frac_up,p_i)

%%

syms p1 p2 p3 p4 p5 a1 a2 a3 a4 a5 a6 a7 a8 a9

A=[a9 0 0 a3 a6;0 0 0 0 0;0 0 0 0 0;a7 0 0 a1 a2; a8 0 0 a4 a5];
bb=[p1,p2,p3,p4,p5]*A*[p1;p2;p3;p4;p5];
collect(bb)

