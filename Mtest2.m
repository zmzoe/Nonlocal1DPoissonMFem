
format short


Value=zeros(4,8);
for l=2:2
delta=0.5;
m=2^l;

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
    
    
    
    %%  ~~~~~~\|G_delta\|^2 ~~~~~~ %%
    
    GD=sparse(n,n);
    
    for i=1:n
        v4=1/5*(X(i+1)^5-X(i)^5)*(A(i,:)'*A(i,:));
        v3=1/4*(X(i+1)^4-X(i)^4)*((A(i,:)'*B(i,:))+B(i,:)'*A(i,:));
        v2=1/3*(X(i+1)^3-X(i)^3)*((B(i,:)'*B(i,:))+C(i,:)'*A(i,:)+A(i,:)'*C(i,:));
        v1=1/2*(X(i+1)^2-X(i)^2)*((B(i,:)'*C(i,:))+C(i,:)'*B(i,:));
        v0=h*(C(i,:)'*C(i,:));
        GD=GD+v0+v1+v2+v3+v4;
    end
    
    %% ~~~~~~\|p\|_0^2 ~~~~~~ %%
    
    
    pl_1=2/3*h*ones(n,1);
    pl_2=h/6*ones(n-1,1);
    
    PL=diag(pl_1);
    PL=PL+diag(pl_2,1)+diag(pl_2,-1);
    PL(n,1)=h/6;
    PL(1,n)=h/6;
    
    %{
    ph_1=1/h*ones(n,1);
    ph_2=-1/h*ones(n-1,1);
    
    PH=diag(ph_1);
    PH=PH+diag(ph_2,1)+diag(ph_2,-1);
    PH(n,1)=-1/h;
    PH(1,n)=-1/h; 
    
    PHH=PH+PL;
    eh1=eig(PHH);
    Value(6,l-1)=min(eh1);
    %}
    
    
    %% ~~~~~\|p\|_delta^2~~~~~%%
    
    PLL=PL+GD;
    ex=eig(PLL);
    Value(3,l-1)=min(ex);
    
    %% ~~~~~\|v\|_0^2~~~~~%%
    
    VLL=h*eye(n);
    ey=eig(VLL);
    Value(4,l-1)=min(ey);
    %% (G_delta p, v)
    
    
    Gv=zeros(n,n);
    
    for i=1:n
        t2=1/3*(X(i+1)^3-X(i)^3)*A(i,:);
        t1=1/2*(X(i+1)^2-X(i)^2)*B(i,:);
        t0=h*C(i,:);

        Gv(i,:)=Gv(i,:)+t0+t1+t2;
    end
    
    sb=svd(Gv);
    %sb(sb<=10^(-10))=[];
    Value(2,l-1)=sb(n-2);
    
    %% singular value
    
   SX=chol(PLL);
   SY=chol(VLL);
   
   E=inv(SY)'*Gv*inv(SX);

   sigma = svd(E);
   sigma(sigma<=10^(-10))=[];
   %s_min=sigma(n-1);
   %s_min=min(sigma);
   
   Value(1,l-1)=sigma(n-2);
end


figure
plot([2^4,2^5,2^6,2^7,2^8,2^9,2^10,2^11].^(-1), Value(1,:),'--s','LineWidth',2);

set(gca, 'LineWidth', 1.5);
legend('MinSinValue');
title('MinSinvalues and h(\delta=0.5)')
xlabel('h')
ylabel('MinSinvalue')
 
figure
plot([2^4,2^5,2^6,2^7,2^8,2^9,2^10,2^11], Value(1,:),'--s','LineWidth',2);

set(gca, 'LineWidth', 1.5);
legend('MinSinValue');
title('MinSinvalues and n(\delta=0.5)')
xlabel('#points n+1')
ylabel('MinSinvalue')

figure
plot([2^4,2^5,2^6,2^7,2^8,2^9,2^10,2^11], Value(1,:),'--s','LineWidth',2);
hold on
plot([2^4,2^5,2^6,2^7,2^8,2^9,2^10,2^11], Value(2,:),'-','LineWidth',2);
hold on
plot([2^4,2^5,2^6,2^7,2^8,2^9,2^10,2^11], Value(3,:),'-o','LineWidth',2);
hold on
plot([2^4,2^5,2^6,2^7,2^8,2^9,2^10,2^11], Value(4,:),'-+','LineWidth',2);

set(gca, 'LineWidth', 1.5);
legend('MinSingularValue','MinSinValueB','MinEigValueXX','MinEigValueYY');
title('MinSinEigvalues and n(\delta=0.5)')
xlabel('#points(1/h)')
ylabel('MinEig-or-Sinvalue')


figure
plot([2^4,2^5,2^6,2^7,2^8,2^9,2^10,2^11].^(-1), Value(1,:),'--s','LineWidth',2);
hold on
plot([2^4,2^5,2^6,2^7,2^8,2^9,2^10,2^11].^(-1), Value(2,:),'-','LineWidth',2);
hold on
plot([2^4,2^5,2^6,2^7,2^8,2^9,2^10,2^11].^(-1), Value(3,:),'-o','LineWidth',2);
hold on
plot([2^4,2^5,2^6,2^7,2^8,2^9,2^10,2^11].^(-1), Value(4,:),'-+','LineWidth',2);

set(gca, 'LineWidth', 1.5);
legend('MinSingularValue','MinSinValueB','MinEigValueXX','MinEigValueYY');
title('MinSinEigvalues and h(\delta=0.5)')
xlabel('h')
ylabel('MinEig-or-Sinvalue')
%{
figure

plot([2^3,2^4,2^5,2^6,2^7,2^8,2^9], Value(6,:),'-+','LineWidth',2);

set(gca, 'LineWidth', 1.5);
legend('MinEigH1');
text(1250,0.6,'/delta=0.5')
text(1250,0.4,'h=/delta/m')
title('Relationship of MinSinvalues and h')
xlabel('#points(1/h)')
ylabel('MinSinvalue')   
%}
%{
figure

plot([2^3,2^4,2^5,2^6], Value(5,:),'-+','LineWidth',2);

set(gca, 'LineWidth', 1.5);
legend('MinSinValueLocalH1');
text(1250,0.6,'/delta=0.5')
text(1250,0.4,'h=/delta/m')
title('Relationship of MinEigvalues and h')
xlabel('#points(1/h)')
ylabel('MinEigvalue')  
%}

