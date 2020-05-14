%   Main
% 
%   Nonlocal formulation with local penalty term
%%
% old definition of G
%    p+Gu = g
%    G'p = f-------Gg
%
%  (p,q)+(Gq,u)+sum [q]lambda=(g,q)
%   (Gp,v)=(f,v)-gamma h^rho
%   sum [p]mu=gamma h^rho
%%

%% Parameter
format long 

delta  =   0.5;
m      =  4;
h      =   delta/m;
n      =   1/h;
gamma  =   0;
rho    =   2;

%% basis 1 discontinuous P1- P0

%% P1


Alpha = eye(2*n, 2*n);


%% P0

Alpha_p0 = zeros(n, 2*n);

for i = 1:n
    Alpha_p0(i,2*i-1)   = 1;
    Alpha_p0(i,2*i)     = 1;
 
end


%% Linear system

A = zeros(2*n,2*n);

for i = 1:2*n
    for j = 1:2*n
         A(i,j) = innerproduct(Alpha(i,:),Alpha(j,:), delta, m );       
    end
end


B = zeros(2*n,n);

for i = 1:2*n
    for j = 1:n
         B(i,j) = Pone(Alpha(i,:),Alpha_p0(j,:), delta, m );       
    end
end



MM=zeros(3*n,3*n);
MM(1:2*n,1:2*n)=A;
MM(1:2*n,2*n+1:3*n)=B;
MM(2*n+1:3*n,1:2*n)=B';

CC=zeros(2*n,n);
CC(1,n)     =   1;
CC(2,1)     =  -1;%
for k  =   2:n
    CC(2*k-1,k-1)   =   1;
    CC(2*k,k)       =  -1;%
end


    Mnew=zeros(4*n,4*n);
    Mnew(1:3*n,1:3*n)=MM;
    Mnew(1:2*n,3*n+1:4*n)=-CC;
    Mnew(3*n+1:4*n,1:2*n)=-CC';





    ureal = @(x) sin(2*pi*x);  
    preal = @(x) 2*pi*cos(2*pi*x);   
    
    f     = @(x) (4*sin(2*pi*x)*(cos(2*pi*delta)/2 - delta^2*pi^2 + delta*pi*sin(2*pi*delta) - 1/2))/(delta^4*pi^2);
    
    Gg    = @(x) (2*(sin(2*pi*x) + cos(2*pi*delta)*sin(2*pi*x) + sin(2*pi*delta)*cos(2*pi*x) + sin(2*pi*x)/(delta^2*pi^2) - 2*delta*pi*cos(2*pi*x) - (cos(2*pi*delta)*sin(2*pi*x))/(delta^2*pi^2) - (2*sin(2*pi*delta)*sin(2*pi*x))/(delta*pi)))/delta^2;
    %按照原来定义算的
 
    g     = @(x) 2*pi*cos(2*pi*x) - (2*((cos(2*pi*x) - cos(2*pi*(delta - x)))/(2*pi) + delta*sin(2*pi*x)))/delta^2;

    HH     = @(x) f(x)+Gg(x);
    
    F     =   Fgauss(h,HH); 
     
    Fhere =   F+gamma*h^rho*ones(n,1) ;
    
    
    G     =   GgaussDirect(h,g);
    FF    =   [G;Fhere];
    
    Fnew  =[FF;gamma*h^rho*ones(n,1)];
    
    
    
    
    
    
   % [xx,flag,relres] = gmres(MM,FF); no penalty
                      
     [xx,flag,relres] = gmres(Mnew,Fnew); 
       %[xx,flag,relres]=gmres(Mnew,Fnew,[],[],[],[],[],x0);
    
    
    pp=xx(1:2*n);
   
   
   
    uu=xx(2*n+1:3*n);
    
     
    errP=getDeltaError_dgpenalty(Alpha,pp,f,Gg,delta,m);
    errP=sqrt(errP);
    
    errP_L2=getL2Error_ppenalty(pp,preal,delta,m);
    errP_L2=sqrt(errP_L2);
    
    errU=getL2Error(uu,ureal,delta,m);
    errU=sqrt(errU);
   
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
    
    
   
    
   
    
    
    