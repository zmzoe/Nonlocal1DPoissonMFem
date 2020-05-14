%function [s_min]=MinSingularValue(delta,m,x, s)

delta=0.5;
m=10;

h=delta/m;
n=1/h;
X=0:h:(n+m)*h;
v=sym(zeros(1,(m+2)));



D=sparse(n,n);%\|p\|_delta^2
F=sparse(n,n);%(Gp,v)
tic
for i=1:n
    if (i+m+1)>n
        J=[i:n,1:(m+i+1-n)];
    else
        J=i:(i+m+1);
    end

    v(1)=(X(i+1)-x)^2/(2*h)-m*(X(i+1)-x);
    v(2)=int(phi_iplus(X,h,i+1,x+s),s,0,X(i+1)-x)+int(phi_iminus(X,h,i+1,x+s),s,X(i+1)-x,X(i+2)-x)-int(phi_iplus(X,h,i+1,x),s,0,delta);
    v(m+1)=int(phi_iplus(X,h,i+m,x+s),s,X(i+m-1)-x,X(i+m)-x)+int(phi_iminus(X,h,i+m,x+s),s,X(i+m)-x,delta);
    v(m+2)=int(phi_iplus(X,h,i+m+1,x+s),s,X(i+m)-x,delta);
    for k=3:m      
        v(k)=h;
    end
 
    v=2/(delta^2)*v;
    
    
    VV=int(v'*v,x,X(i),X(i+1));
    
    VV=double(VV);
    
    DD=zeros(n,n);
 
    for p=J
        pp=find(J==p);
        for q=J     
            qq=find(J==q);
            DD(p,q)=DD(p,q)+VV(pp,qq);
        end      
    end

    D=DD+D;   
    
  UU=int(v,x,X(i),X(i+1));
  FF=zeros(n,n);
  for p=J
        pp=find(J==p);
        FF(i,p)=FF(i,p)+UU(pp);       
  end
  F=FF+F;
   
end
toc            
 
tic
AA=zeros(n,n);

  for i=1:(n-1)
      a_i=ai(X,h,i,x);
      for j=i:i+1
          for k=i:i+1
               AA(j,k)=AA(j,k)+a_i(j-i+1,k-i+1);
          end
      end
  end          
     
  a_n=ai(X,h,n,x);  
  
  AA(n,n)=AA(n,n)+a_n(1,1);
  AA(n,1)=AA(n,1)+a_n(2,1);
  AA(1,n)=AA(1,n)+a_n(1,2);
  AA(1,1)=AA(1,1)+a_n(2,2);    
  
  D=D+AA;
  toc

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
W(W<=(10^-8))=[];

s_min=min(W);
