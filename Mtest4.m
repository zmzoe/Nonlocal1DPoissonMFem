


function [pS,Sigma]=Mtest4(delta, n)

    
h=1/n;

X=0:h:(n+2)*h;



A=sparse(n,n);
B=sparse(n,n);
C=sparse(n,n);

va=zeros(1,3);
vb=zeros(1,3);
vc=zeros(1,3);

for i=1:n
    I=i*ones(1,3);
    if i<=(n-2)
        J=(i:i+2);
    else
        J=[i:n,1:2+i-n];
    end
    
    va(1)=1/(2*h);
    vb(1)=-X(i+1)/h+delta/h;
    vc(1)=X(i+1)^2/(2*h)-delta*X(i+1)/h;
    
    
    va(2)=-1/h;
    vb(2)=-(h+delta)/h-delta/h;
    vc(2)=-(X(i+1)*(2*X(i)-X(i+1)))/(2*h)-(delta-X(i+1))*(delta+X(i+1)-2*X(i+2))/(2*h)+delta*X(i)/h;
    
    va(3)=1/(2*h);
    vb(3)=(2*delta - 2*X(i+1))/(2*h);
    vc(3)=(delta - X(i+1))^2/(2*h);
    
    A=A+sparse(I,J,va,n,n);
    B=B+sparse(I,J,vb,n,n);
    C=C+sparse(I,J,vc,n,n);
          
    
end

    A=2/(delta^2)*A;
    B=2/(delta^2)*B;
    C=2/(delta^2)*C;
    
    
    
     VLL=h*eye(n);
    for i=1:n
        VLL(1:n-1,i)=VLL(1:n-1,i)-VLL(n,i);     
    end
    VLL(n,:)=[];
   
    for i=1:n-1
        VLL(i,1:n-1)=VLL(i,1:n-1)-VLL(i,n);     
    end
    VLL(:,n)=[];
    
    
    %% (G_delta p, v)
    
    
    Gv=zeros(n,n);
    
    for i=1:n
        t2=1/3*(X(i+1)^3-X(i)^3)*A(i,:);
        t1=1/2*(X(i+1)^2-X(i)^2)*B(i,:);
        t0=h*C(i,:);

        Gv(i,:)=Gv(i,:)+t0+t1+t2;
    end
    
    for i=1:n
        Gv(1:n-1,i)=Gv(1:n-1,i)-Gv(n,i);     
    end
    Gv(n,:)=[];
   
    for i=1:n-1
        Gv(i,1:n-1)=Gv(i,1:n-1)-Gv(i,n);     
    end
    Gv(:,n)=[];
    
    %% ~~~~~~~~\|Gp\|^2 ~~~~~~~~~~~
    
    
    
    GD=sparse(n,n);
    
    for i=1:n
        v4=1/5*(X(i+1)^5-X(i)^5)*(A(i,:)'*A(i,:));
        v3=1/4*(X(i+1)^4-X(i)^4)*((A(i,:)'*B(i,:))+B(i,:)'*A(i,:));
        v2=1/3*(X(i+1)^3-X(i)^3)*((B(i,:)'*B(i,:))+C(i,:)'*A(i,:)+A(i,:)'*C(i,:));
        v1=1/2*(X(i+1)^2-X(i)^2)*((B(i,:)'*C(i,:))+C(i,:)'*B(i,:));
        v0=h*(C(i,:)'*C(i,:));
        GD=GD+v0+v1+v2+v3+v4;
    end
    
    

    for i=1:n
        GD(1:n-1,i)=GD(1:n-1,i)-GD(n,i);     
    end
    GD(n,:)=[];
   
    for i=1:n-1
        GD(i,1:n-1)=GD(i,1:n-1)-GD(i,n);     
    end
    GD(:,n)=[];
    
   
   
   SX=chol(GD);
   SY=chol(VLL);
   
   E=inv(SY)'*Gv*inv(SX);

   sigma = svd(E);
   sigma(sigma<=10^(-10))=[];
   
   Sigma=min(sigma);
    
   
   
   %%
   
    pl_1=2/3*h*ones(n,1);
    pl_2=h/6*ones(n-1,1);
    
    PL=diag(pl_1);
    PL=PL+diag(pl_2,1)+diag(pl_2,-1);
    PL(n,1)=h/6;
    PL(1,n)=h/6;
    
    for i=1:n
        PL(1:n-1,i)=PL(1:n-1,i)-PL(n,i);     
    end
    PL(n,:)=[];
   
    for i=1:n-1
        PL(i,1:n-1)=PL(i,1:n-1)-PL(i,n);     
    end
    PL(:,n)=[];
    
   
    
    PLL=PL+GD;
    SXX=chol(PLL);
    EE=inv(SY)'*Gv*inv(SXX);

   pS = svd(EE);
   pS(pS<=10^(-10))=[];
   pS=min(pS);
   
%end
    




