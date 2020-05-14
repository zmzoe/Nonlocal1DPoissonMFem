%   Main
% 
%   Nonlocal formulation with local penalty term
%%
% old definition of G
%
%  (p,q)+(Gq,u)+sum [q]lambda=(g,q)
%   (Gp,v)=(f-Gg,v):=(HHv)
%   sum [p]mu=gamma h^rho
%%---------
% G=-G
%%-----------
format long
%% Parameter

delta  =   0.5;
m      =   5;
h      =   delta/m;
n      =   1/h;
gamma  =  0;
rho    =   2;

%% basis 1 discontinuous P1- P0

%% P1
Alpha = zeros(2*n-1, 2*n);

for i=1:2*n-1
    Alpha(i,i)   =  1;
    Alpha(i,i+1) = -1;
end

%% P0

Alpha_p0 = zeros(n-1, 2*n);

for i = 1:n-1
    Alpha_p0(i,2*i-1)     =  1;
    Alpha_p0(i,2*i)       =  1;
    Alpha_p0(i,2*i+1)     = -1;
    Alpha_p0(i,2*i+2)     = -1;
 
end


%% Linear system

A = zeros(2*n-1,2*n-1);

for i = 1:2*n-1
    for j = 1:2*n-1
         A(i,j) = innerproduct(Alpha(i,:),Alpha(j,:), delta, m );       
    end
end


B = zeros(2*n-1,n-1);

for i = 1:2*n-1
    for j = 1:n-1
         B(i,j) = Pone(Alpha(i,:),Alpha_p0(j,:), delta, m );       
    end
end

MM=zeros(3*n-2,3*n-2);
MM(1:2*n-1,1:2*n-1)=A;
MM(1:2*n-1,2*n:3*n-2)=B;
MM(2*n:3*n-2,1:2*n-1)=B';


CC=zeros(2*n,n);
CC(1,n)     =   1;
CC(2,1)     =  -1;
for k  =   2:n
    CC(2*k-1,k-1)   =   1;
    CC(2*k,k)       =  -1;
end

DD=zeros(2*n-1,n);
for i=1:2*n-1
    for j=1:n
        DD(i,j)=Alpha(i,:)*CC(:,j);
    end
    
end


    Mnew=zeros(4*n-2,4*n-2);
    Mnew(1:3*n-2,1:3*n-2)=MM;
    Mnew(1:2*n-1,3*n-1:4*n-2)=-DD;
    Mnew(3*n-1:4*n-2,1:2*n-1)=-DD';





    ureal = @(x) sin(2*pi*x);  
    preal = @(x) 2*pi*cos(2*pi*x);   
    
    f     = @(x) (4*sin(2*pi*x)*(cos(2*pi*delta)/2 - delta^2*pi^2 + delta*pi*sin(2*pi*delta) - 1/2))/(delta^4*pi^2);
    
    Gg    = @(x) (2*(sin(2*pi*x) + cos(2*pi*delta)*sin(2*pi*x) + sin(2*pi*delta)*cos(2*pi*x) + sin(2*pi*x)/(delta^2*pi^2) - 2*delta*pi*cos(2*pi*x) - (cos(2*pi*delta)*sin(2*pi*x))/(delta^2*pi^2) - (2*sin(2*pi*delta)*sin(2*pi*x))/(delta*pi)))/delta^2;
    %按照原来定义算的
 
    g     = @(x) 2*pi*cos(2*pi*x) - (2*((cos(2*pi*x) - cos(2*pi*(delta - x)))/(2*pi) + delta*sin(2*pi*x)))/delta^2;

    HH     = @(x) f(x)+Gg(x);
    
    gq  =  zeros(2*n-1,1);
    for j=1:2*n-1
        gq(j,1)=Gauss_anybase(Alpha(j,:),g,h); %(g,q)
    end
    
    HHv  =  zeros(n-1,1);
    for j=1:n-1
        HHv(j,1)=Gauss_anybase(Alpha_p0(j,:),HH,h); %(HH,v)
    end
    
  %  F     =   [gq;HHv;zeros(n,1)]; 
   
    
    
    F     =   [gq;HHv;gamma*h^rho*ones(n,1)]; 
    
    
  
    
  
    % [xx,flag,relres] = gmres(Mnew,F); 
     xx=Mnew\F;
    
    
    pp1=xx(2:2*n-1)-xx(1:2*n-2);
    pp=[xx(1);pp1;-xx(2*n-1)];
   
    ppn=[pp;pp(1)];
    
   
    uu0=xx(2*n:3*n-2);
    uu1=uu0(2:n-1)-uu0(1:n-2);
    uu=[uu0(1);uu1;-uu0(n-1)];
    
     
    
    
     errP=getDeltaError_anybase_dg(Alpha,xx,HH,delta,m);
     errP=sqrt(errP);%修改
     
    
     
     
     
     err_P_L2=getL2Error_local_p_dg(pp,preal,h);
     err_P_L2=sqrt(err_P_L2);
    
     err_U=getL2Errorlocal(uu,ureal,h);
     err_U=sqrt(err_U);
    
    
    
    
    
    
    
   %% Figure 1------- P ---------P_h---------- P-P_h
    
    
    
    figure 
   for i=1:n
        sx=[(i-1)*h,i*h];
        sy=[pp(2*i-1),pp(2*i)];
        sz=[preal(sx(1)),preal(sx(2))];
        
        o1=line(sx,sz,'Color',[1 1 0],'LineWidth',4);
        hold on
        o2=line(sx,sy,'Color',[0 0 0],'LineWidth',1);
        hold on
        o3=plot((sx(1)+sx(2))/2,(sy(1)+sy(2))/2,'k*');
        hold on
   end     
    legend([o1,o2,o3],{ 'realP','numericalP','numPmidpoint'})
    xlabel('x')
    ylabel('function value')
    title('NumSolP')
    
    
    
   
    
    figure 
    for i=1:n
        sx=[(i-1)*h,i*h];
        sy=[uu(i),uu(i)];
        sz=[ureal(sx(1)),ureal(sx(2))];
        o1=line(sx,sz,'Color',[1 1 0],'LineWidth',4);
        hold on
        o2=line(sx,sy,'Color',[0 0 0],'LineWidth',1);
        hold on
        o3=plot((sx(1)+sx(2))/2,(sy(1)+sy(2))/2,'k*');
        hold on
    end
    legend([o1,o2,o3],{'realU','numericalU','midpoint'})
    xlabel('x')
    ylabel('function value')
    title('NumSolU')
    
    
   
    
   
    
    
    