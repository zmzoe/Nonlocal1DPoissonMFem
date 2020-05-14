
clear
clc
syms x x_1 x_n x_n x_n1 s h t delta k n m i x_3 x_4 x_n2 x_2 x_3 x_n3 x_i x_i1 x_i2 x_m2 x_m1 x_i_1


%{
ff1=(x+s-(n+t-2)*h)/h;
ff2=((n+t)*h-x-s)/h;

f1=(x-(t-2)*h)/h;
f2=(t*h-x)/h;

%}
%


f_iplus=(x-(t-2)*h)/h;
f_i1plus=(x-(t-1)*h)/h;
f_iminus=(t*h-x)/h;
f_i1minus=((t+1)*h-x)/h;


F=int(f_iplus^2,x,(t-2)*h,(t-1)*h)+int(f_iminus^2,x,(t-1)*h,(t)*h);

FIJ=int(f_iminus*f_i1plus,x,(t-1)*h,(t)*h);



