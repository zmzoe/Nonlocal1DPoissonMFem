
function [alpha_i,beta_i,gamma_i]=gradp(delta,m,II)

h=delta/m;
n=1/h;
X=0:h:(n+m)*h;




A=sparse(n,n);%G_delta
B=sparse(n,n);
C=sparse(n,n);

va=zeros(1,m+2);
vb=zeros(1,m+2);
vc=zeros(1,m+2);

for i=1:n
    I=i*ones(1,m+2);
   
    
    if (i+m+1)<=n
        J=i:(i+m+1);
    else
        J=[i:n,1:(m+i+1-n)];
    end
    
    if m==1
        va(1)=1/(2*h);
        vb(1)=delta/h - X(i+1)/h;
        vc(1)=X(i+1)^2/(2*h) - (delta*X(i+1))/h;
        
        va(2)=-1/h;
        vb(2)= (2*X(i) - X(i+1))/(2*h) - (delta + X(i+1) - 2*X(i+2))/(2*h) - delta/h + X(i+1)/(2*h) - (delta - X(i+1))/(2*h);
        %vb(2)=(X(i)+X(i+2))/h-2;
        vc(2)=(delta*X(i))/h - (X(i+1)*(2*X(i) - X(i+1)))/(2*h) - ((delta - X(i+1))*(delta + X(i+1) - 2*X(i+2)))/(2*h);
        %vc(2)=-(X(i+1)*(2*X(i) - X(i+1)))/(2*h) - ((h - X(i+1))*(h + X(i+1) - 2*X(i+2)))/(2*h)+X(i);
        
        va(3)=1/(2*h);
        vb(3)=(2*delta - 2*X(i+1))/(2*h);
        vc(3)=(delta - X(i+1))^2/(2*h);
        
    else
    
      va(1)=1/(2*h);
      vb(1)=-1/h*X(i+1)+delta/h;
      vc(1)=1/(2*h)*X(i+1)^2-delta*X(i+1)/h;
    
      va(2)=-1/(2*h);
      vb(2)=1/h*X(i)-delta/h;
      vc(2)=1/(2*h)*X(i+1)^2-X(i)*X(i+1)/h+h/2+delta*X(i)/h;
    
      va(m+1)=-1/(2*h);
      vb(m+1)=X(i+m+1)/h-X(i+m)/h-(delta-X(i+m))/h;
      vc(m+1)=h/2-(delta^2-X(i+m)^2)/(2*h)+(X(i+m+1)*(delta-X(i+m))/h);
    
    
      va(m+2)=1/(2*h);
      vb(m+2)=-1/h*X(i+m)+delta/h;
      vc(m+2)=1/(2*h)*(delta^2-X(i+m)^2)-X(i+m)*(delta-X(i+m))/h;
    
      vc(3:m)=h;
    end                   
    
    A=A+sparse(I,J,va,n,n);
    B=B+sparse(I,J,vb,n,n);
    C=C+sparse(I,J,vc,n,n); 
    
end

%% ~~~~~~ G_delta^* p over (x_i,x_i1) ~~~~~~%%

    
    
    A=2/(delta^2)*A;
    B=2/(delta^2)*B;
    C=2/(delta^2)*C;
    
    alpha_i = A(II,:);
    beta_i  = B(II,:);
    gamma_i = C(II,:);
    
   

end
