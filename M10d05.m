%h=delta/M
%Omega=(0,1) periodical boundary condition

%delta=0.5
%M=10
%h=0.05
%#points=21

%main

%%
clc;
format long

delta=0.5;
m=10;
h=delta/m;
n=1/h;
X=0:h:(n+m)*h;


%%
syms p_1 p_2 p_3 p_4 p_5 p_6 p_7 p_8 p_9 p_10 p_11 p_12 p_13 p_14 p_15 p_16 p_17 p_18 p_19 p_20 s x
syms v_1 v_2 v_3 v_4 v_5 v_6 v_7 v_8 v_9 v_10 v_11 v_12 v_13 v_14 v_15 v_16 v_17 v_18 v_19 v_20

P=sym(zeros(n+1,1));

for i=1:n
    P(i)=p_i;
end
P=[p_1 p_2 p_3 p_4 p_5 p_6 p_7 p_8 p_9 p_10 p_11 p_12 p_13 p_14 p_15 p_16 p_17 p_18 p_19 p_20 p_1];

P_1=[p_1 p_2 p_3 p_4 p_5 p_6 p_7 p_8 p_9 p_10 p_11 p_12 p_13 p_14 p_15 p_16 p_17 p_18 p_19 p_20];
V=[v_1 v_2 v_3 v_4 v_5 v_6 v_7 v_8 v_9 v_10 v_11 v_12 v_13 v_14 v_15 v_16 v_17 v_18 v_19 v_20];

%% ~~~~~~~~~~xi < x < xi+1~~~~~~~~~~~~~~~G_delta^* . d=delta  .             

I11=0;
I22=0;
I33=0;
for i=1:n %x_i
    
    g_ii=g(P,i,X,h,x);
    f_ii=int(g_ii,s,0,delta); % \int_0^delta p(x) ds
    
    g_i1=g(P,i,X,h,x+s); %1st -term
    g_im=g(P,i+m,X,h,x+s); %m+1 th - term
    
    f_jj=int(g_im,s,X(i+m)-x,delta)+int(g_i1,s,0,X(i+1)-x);
    
    for j=1:m-1
    g_i=g(P,i+j,X,h,x+s);
    f_j=int(g_i,s,X(i+j)-x,X(i+j+1)-x);
    f_jj=f_j+f_jj;
    end
    f_jj=f_jj-f_ii;
    f_jj=2/(delta^2)*f_jj;
    %G_delta finished
    
    I11=I11+int(f_jj^2,x,X(i),X(i+1)); %\|G_delta p\|^2
    I22=I22+int(g_ii^2,x,X(i),X(i+1)); %\|p\|^2
    
    I33=I33+int(f_jj,x,X(i),X(i+1))*V(i); % (G_delta p, v)
end
    
     I1=I11+I22; %\|p\|_delta^2
     
    % collect(I1,p_1)
     
     I2=expand(I1); %simplify


     M=1/2*hessian(I2); %this is the M  p*M*p^T=\|p\|_delta^2

     
%% ~~~~~~~~~~~~~~~ \|v\|_0^2 ~~~~~~~~~~~~~~~~~~   
     
N=h*eye(n);
     
%% ~~~~~~~~~~~~~~~~~SX & SY ~~~~~~~~~~~~~~~~~~~

SX=chol(M);
SY=chol(N);

%% ~~~~~~~~~~How to compute B~~~~~~~~~~~~~~~~~~~~

I33=expand(I33);

 for i=1:n
     for j=1:n
         B(i,j)=diff(diff(I33,P_1(i)),V(j));
     end
 end
    
TT=inv(SY)'*B'*inv(SX);
%W=svds(TT);
%s_min=min(s);
%W=s_min;[U,S,V] = svd(A)

W = svd(TT);

W(W==0)=[];

s_min=min(W);

