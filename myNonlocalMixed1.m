

function [err_P_delta,err_P_L2,err_U,ERR_u0,ERR0,ERR1]=myNonlocalMixed1(delta,m)

% Parameter
h      =   delta/m;
n      =   1/h;





% basis 1 continuous P1- P0
% P1
Alpha = zeros(n, 2*n);
for i = 2:n
    Alpha(i,2*(i-1)) =     1;
    Alpha(i,2*i-1)   =     1;
end
Alpha(1,1)   =      1;
Alpha(1,2*n) =      1;
% P0
Alpha_p0 = zeros(n, 2*n);
for i = 1:n
    Alpha_p0(i,2*i-1)     =  1;
    Alpha_p0(i,2*i)       =  1;
 
end


% Linear system
A = zeros(n,n);
for i = 1:n
    for j = 1:n
         A(i,j) = innerproduct(Alpha(i,:),Alpha(j,:), delta, m );       
    end
end
B = zeros(n,n);
for i = 1:n
    for j = 1:n
         B(i,j) = Pone(Alpha(i,:),Alpha_p0(j,:), delta, m );       
    end
end
MM=zeros(2*n,2*n);
MM(1:n,1:n)=A;
MM(1:n,n+1:2*n)=-B;
MM(n+1:2*n,1:n)=-B';




% Real solution and RHS
ureal = @(x) sin(6*pi*x); 
preal = @(x) -(2*((cos(6*pi*x) - cos(6*pi*(delta - x)))/(6*pi) + delta*sin(6*pi*x)))/delta^2;
f     = @(x) -(2*(2*sin(6*pi*x) - ((cos(6*pi*delta).*sin(6*pi*x))/9 - sin(6*pi*x)/9 + (2*delta*pi*sin(6*pi*delta).*sin(6*pi*x))/3)/(delta^2*pi^2)))/delta^2;
HHv   =  zeros(n,1);
for j=1:n
    HHv(j,1)=Gauss_anybase(Alpha_p0(j,:),f,h); %(f,v)
end  
F     =   [zeros(n,1);HHv]; 
  



% Gmres solve the system 
[xx,~,~] = gmres(MM,F);  
p_h=[xx(1:n);xx(1)];
u_h=xx(n+1:2*n);

p_r=preal(linspace(0,1-h,n)');
u_r=ureal(linspace(0,1-h,n)');


% Stardard Error 
 err_P_delta=getDeltaError4anybase(Alpha,p_h,f,delta,m);
 err_P_delta=sqrt(err_P_delta);

err_P_L2=getL2Error4p_ph(p_h,preal,h);
err_P_L2=sqrt(err_P_L2);

err_U=getL2Errorlocal(u_h,ureal,h);
err_U=sqrt(err_U);
      



% Superconvergence Error 

coe_u=u_r-u_h(1:n);
ERR_u0=sqrt(coe_u'*coe_u);

coe=p_r-p_h(1:n);
ERR0=sqrt(coe'*A*coe);

ERR1=get_Ip_ph(Alpha,p_h,p_r,delta,m);
ERR1=sqrt(ERR1);





% figure
figure
plot(0:h:1,p_h,'LineWidth',2);
hold on
plot(0:h:1,preal(0:h:1),'LineWidth',2);
    
legend('numerical-p', 'real-p')
xlabel('x')
ylabel('function value')
title('NumSol-p')
    
figure 
for i=1:n
    sx=[(i-1)*h,i*h];
    sy=[u_h(i),u_h(i)];
    sw=[ureal(sx(1)),ureal(sx(2))]-sy;
    o1=line(sx,sy,'LineWidth',2);
    hold on
    o3=line(sx,sw,'Color',[1 0.8 1],'LineWidth',2);
    hold on
end
    hold on
    o2=plot(0:h:1,ureal(0:h:1),'Color',[1 0.5 0],'LineWidth',1);
    legend([o1,o2,o3],{'numerical-u', 'real-u','errU'})
    xlabel('x')
    ylabel('function value')
    title('NumSol-U')
    

  

    
    
    
    
    
    