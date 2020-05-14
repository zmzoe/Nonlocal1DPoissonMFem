format short


Value=zeros(5,6);
for l=3:3
delta=0.2;
m=4;

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
    
    A=A+sparse(I,J,va,n,n);
    B=B+sparse(I,J,vb,n,n);
    C=C+sparse(I,J,vc,n,n); 
    
  
   
end

%% ~~~~~~ G_delta^* p over (x_i,x_i1) ~~~~~~%%

    
    
    A=2/(delta^2)*A;
    B=2/(delta^2)*B;
    C=2/(delta^2)*C;
   
    
    

    %% ~~~~~\|v\|_0^2~~~~~%%
    
    VLL=h*eye(n);
    for i=1:n
        VLL(1:n-1,i)=VLL(1:n-1,i)-VLL(n,i);     
    end
    VLL(n,:)=[];
   
    for i=1:n-1
        VLL(i,1:n-1)=VLL(i,1:n-1)-VLL(i,n);     
    end
    VLL(:,n)=[];
    
    ey=eig(VLL);
    Value(4,l-2)=min(ey);
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
    
    sb=svd(Gv);
    sb(sb<=10^(-10))=[];
    Value(2,l-2)=min(sb);
    
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
    
    ex=eig(full(GD));
    Value(3,l-2)=min(ex);
   
   SX=chol(GD);
   SY=chol(VLL);
   
   E=inv(SY)'*Gv*inv(SX);

   sigma = svd(E);
   sigma(sigma<=10^(-10))=[];
   
   
   Value(1,l-2)=min(sigma); 
   
   
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
    
    ep=eig(full(PL));
    Value(5,l-2)=min(ep);
    
end
   

%{
figure
plot([2^4,2^5,2^6,2^7,2^8,2^9].^(-1), Value(5,:),'--s','LineWidth',2);

set(gca, 'LineWidth', 1.5);
legend('MinEigValueOFpl2');
title('MinSinvalues and h(\delta=0.5)')
xlabel('h')
ylabel('MinSinvalue')


figure
plot([2^4,2^5,2^6,2^7,2^8,2^9], Value(1,:),'--s','LineWidth',2);

set(gca, 'LineWidth', 1.5);
legend('MinSinValue');
title('MinSinvalues and n(\delta=0.5)')
xlabel('#points n+1')
ylabel('MinSinvalue')

figure
plot([2^4,2^5,2^6,2^7,2^8,2^9], Value(1,:),'--s','LineWidth',2);
hold on
plot([2^4,2^5,2^6,2^7,2^8,2^9], Value(2,:),'-','LineWidth',2);
hold on
plot([2^4,2^5,2^6,2^7,2^8,2^9], Value(3,:),'-o','LineWidth',2);
hold on
plot([2^4,2^5,2^6,2^7,2^8,2^9], Value(4,:),'-+','LineWidth',2);

set(gca, 'LineWidth', 1.5);
legend('MinSingularValue','MinSinValueB','MinEigValueXX','MinEigValueYY');
title('MinSinEigvalues and n(\delta=0.5)')
xlabel('#points(1/h)')
ylabel('MinEig-or-Sinvalue')


figure
plot([2^4,2^5,2^6,2^7,2^8,2^9].^(-1), Value(1,:),'--s','LineWidth',2);
hold on
plot([2^4,2^5,2^6,2^7,2^8,2^9].^(-1), Value(2,:),'-','LineWidth',2);
hold on
plot([2^4,2^5,2^6,2^7,2^8,2^9].^(-1), Value(3,:),'-o','LineWidth',2);
hold on
plot([2^4,2^5,2^6,2^7,2^8,2^9].^(-1), Value(4,:),'-+','LineWidth',2);

set(gca, 'LineWidth', 1.5);
legend('MinSingularValue','MinSinValueB','MinEigValueXX','MinEigValueYY');
title('MinSinEigvalues and h(\delta=0.5)')
xlabel('h')
ylabel('MinEig-or-Sinvalue')    
    
%} 