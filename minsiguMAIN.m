%h=delta/M
%Omega=(0,1) periodical boundary condition

%delta=0.5
%M=10
%h=0.05
%#points=21

%main

%%
clear;
clc;

format short

delta=0.5;
m=10;
h=delta/m;
n=1/h;
X=0:h:(n+m)*h;



syms x s

%% ~~~~~~~~~~xi < x < xi+1~~~~~~~~~~~~~~~G_delta^* p  

A=sym(zeros(n,n)); %G_delta(x_i, x_i+1)
D=sym(zeros(n,n)); %int ((G)^2)
F=sym(zeros(n,n)); %int(G_delta)


tic
for i=1:n
    
    A(i,i)=int(phi_iminus(X,h,i,x+s),s,0,X(i+1)-x)-int(phi_iminus(X,h,i,x),s,0,delta);
    
    if i<n
        A(i,i+1)=int(phi_iplus(X,h,i+1,x+s),s,0,X(i+1)-x)+int(phi_iminus(X,h,i+1,x+s),s,X(i+1)-x,X(i+2)-x)-int(phi_iplus(X,h,i+1,x),s,0,delta);
    else
        A(i,i+1-n)=A(i,i+1-n)+int(phi_iplus(X,h,i+1,x+s),s,0,X(i+1)-x)+int(phi_iminus(X,h,i+1,x+s),s,X(i+1)-x,X(i+2)-x)-int(phi_iplus(X,h,i+1,x),s,0,delta);
    end
    
    if i<(n-m+1)
        A(i,i+m)=int(phi_iplus(X,h,i+m,x+s),s,X(i+m-1)-x,X(i+m)-x)+int(phi_iminus(X,h,i+m,x+s),s,X(i+m)-x,delta);
    else
        A(i,i+m-n)=A(i,i+m-n)+int(phi_iplus(X,h,i+m,x+s),s,X(i+m-1)-x,X(i+m)-x)+int(phi_iminus(X,h,i+m,x+s),s,X(i+m)-x,delta);
    end
    
    if i<(n-m)
        A(i,i+m+1)=int(phi_iplus(X,h,i+m+1,x+s),s,X(i+m)-x,delta);
    else
        A(i,i+m+1-n)=A(i,i+m+1-n)+int(phi_iplus(X,h,i+m+1,x+s),s,X(i+m)-x,delta);
    end
    
    for j=(i+2):(i+m-1)
        
        if j< n
             A(i,j)=int(phi_iplus(X,h,j,x+s),s,X(j-1)-x,X(j)-x)+int(phi_iminus(X,h,j,x+s),s,X(j)-x,X(j+1)-x);
         else
             A(i,j-n+1)=A(i,j-n+1)+int(phi_iplus(X,h,j,x+s),s,X(j-1)-x,X(j)-x)+int(phi_iminus(X,h,j,x+s),s,X(j)-x,X(j+1)-x);    
        end   
    end
end

toc

    A=2/(delta^2)*A;
    
  %% ~~~~~~~~~~~\|p\|^2
  
  AA=zeros(n,n);
  tic
  for i=1:(n-1)
      a_i=ai(X,h,i,x);
      for j=i:i+1
          for k=i:i+1
               AA(j,k)=AA(j,k)+a_i(j-i+1,k-i+1);
          end
      end
  end
  toc
  
  a_n=ai(X,h,n,x);  
  
  AA(n,n)=AA(n,n)+a_n(1,1);
  AA(n,1)=AA(n,1)+a_n(2,1);
  AA(1,n)=AA(1,n)+a_n(1,2);
  AA(1,1)=AA(1,1)+a_n(2,2);
  
  
  
                      
              
  %% ~~~~~~~~~~~~~~~~~~~~~~~~D~~~~~~~~~~~~~~~~~~~~~~
  
  %{
    MM=sym(zeros(n,n));
   tic 
    for i=1:n
        M=A(i,:)'*A(i,:);
        for j=1:n
            for k=1:n
                MM(j,k)=int(M(j,k),x,X(i),X(i+1));
            end
        end
        D=D+MM;
    end
    toc
    %}
  
  
  for i=1:n
        M=A(i,:)'*A(i,:);
        M=int(M,x,X(i),X(i+1));
        D=D+M;
  end
%% ~~~~~~~~~~~`\|p\|_delta^2~~~~~~~~~~~~~~~~~




D=D+AA;




D=double(D);
%% ~~~~~~~~~~~~~~~~~~~~ F ~~~~~~~~~~~~~~~~~~

FF=sym(zeros(n,n));

for i=1:n
    for j=1:n
        FF(i,j)=int(A(i,j),x,X(i),X(i+1));      
    end
    F=F+FF;
end

F=double(F);

%%

%% ~~~~~~~~~~~~~~~ \|v\|_0^2 ~~~~~~~~~~~~~~~~~~   
     
N=h*eye(n);
     
%% ~~~~~~~~~~~~~~~~~SX & SY ~~~~~~~~~~~~~~~~~~~
tic
SX=chol(D);
SY=chol(N);
toc


    
TT=inv(SY)'*F*inv(SX);

tic
W = svd(TT);
toc
W(W<=(10^-14))=[];

s_min=min(W);




