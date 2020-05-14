



function [alpha_i,beta_i,gamma_i]=gradp_dg(delta,m,II)

h=delta/m;
n=1/h;
X=0:h:(n+m)*h;




A_odd=sparse(n,2*n);%G_delta
B_odd=sparse(n,2*n);
C_odd=sparse(n,2*n);

A_even=sparse(n,2*n);%G_delta
B_even=sparse(n,2*n);
C_even=sparse(n,2*n);

va_odd=zeros(1,m+1);
vb_odd=zeros(1,m+1);
vc_odd=zeros(1,m+1);

va_even=zeros(1,m+1);
vb_even=zeros(1,m+1);
vc_even=zeros(1,m+1);

JJ_odd  = zeros(1,m+1);
JJ_even = zeros(1,m+1);

for i=1:n
    I=i*ones(1,m+1);
   
    
    if i+m <= n
        J = i : i+m;
    else
        J=[i:n,1:(i+m-n)];
    end
   
    
    for j = 1: m+1
        JJ_odd(j)  = 2*J(j)-1;
        JJ_even(j) = 2*J(j);
    end
    
    
      va_odd(1)=1/(2*h);
      vb_odd(1)=-1/h*X(i+1)+delta/h;
      vc_odd(1)=1/(2*h)*X(i+1)^2-delta*X(i+1)/h;
      
      va_even(1)=-1/(2*h);
      vb_even(1)=1/h*X(i)-delta/h;
      vc_even(1)=delta*X(i)/h-X(i+1)*(2*X(i)-X(i+1))/(2*h);
    
    
      va_odd(m+1)=-1/(2*h);
      vb_odd(m+1)=-(delta+X(i+m)-2*X(i+m+1))/(2*h)-(delta-X(i+m))/(2*h);
      vc_odd(m+1)=-(delta-X(i+m))*(delta+X(i+m)-2*X(i+m+1))/(2*h);
    
      
      va_even(m+1)=1/(2*h);
      vb_even(m+1)=2*(delta-X(i+m))/(2*h);
      vc_even(m+1)=(delta-X(i+m))^2/(2*h);
    
      
    
      vc_odd(2:m)=h/2;
      vc_even(2:m)=h/2;
   
    A_odd=A_odd+sparse(I,JJ_odd,va_odd,n,2*n);
    B_odd=B_odd+sparse(I,JJ_odd,vb_odd,n,2*n);
    C_odd=C_odd+sparse(I,JJ_odd,vc_odd,n,2*n); 
    
    A_even=A_even+sparse(I,JJ_even,va_even,n,2*n);
    B_even=B_even+sparse(I,JJ_even,vb_even,n,2*n);
    C_even=C_even+sparse(I,JJ_even,vc_even,n,2*n); 
    
end

%% ~~~~~~ G_delta^* p over (x_i,x_i1) ~~~~~~%%

    A=A_odd+A_even;
    B=B_odd+B_even;
    C=C_odd+C_even;
    
    A=2/(delta^2)*A;
    B=2/(delta^2)*B;
    C=2/(delta^2)*C;
    
    alpha_i = A(II,:);
    beta_i  = B(II,:);
    gamma_i = C(II,:);
    
   

end
