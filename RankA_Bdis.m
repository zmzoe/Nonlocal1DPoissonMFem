% Main
% Rank discontinuous A B
% p u and f is derived from the nonlocal operator
%% Parameter

delta  =   0.5;
m      =   4;
h      =   delta/m;
n      =   1/h;


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

RA=rank(A);
SA=size(A);
RB=rank(B);
SB=size(B);

MM=zeros(3*n-2,3*n-2);
MM(1:2*n-1,1:2*n-1)=A;
MM(1:2*n-1,2*n:3*n-2)=-B;
MM(2*n:3*n-2,1:2*n-1)=-B';


    ureal = @(x) sin(2*pi*x);  
    preal = @(x) -(2*((cos(2*pi*x) - cos(2*pi*(delta - x)))/(2*pi) + delta*sin(2*pi*x)))/delta^2;    
    f     = @(x) (4*sin(2*pi*x)*(cos(2*pi*delta)/2 - delta^2*pi^2 + delta*pi*sin(2*pi*delta) - 1/2))/(delta^4*pi^2);

    F  =   Fgauss(h,f); 
    
    Fhere = F(1:n-1,1)-F(2:n,1);
    
    FF    = [zeros(2*n-1,1);Fhere];
    
    
    xx=MM\FF;
    %[xx,flag,relres] = minres(MM,FF);
    
    pp1=xx(2:2*n-1)-xx(1:2*n-2);
    pp=[xx(1);pp1;-xx(2*n-1)];
   
    
   
    uu0=xx(2*n:3*n-2);
    uu1=uu0(2:n-1)-uu0(1:n-2);
    uu=[uu0(1);uu1;-uu0(n-1)];
   
    
    err=getError_sin_ph(Alpha,xx,preal,delta,m);
    
    
    errP=getDeltaError_dgrank(Alpha,xx,f,delta,m);
    errP=sqrt(errP);
    
    
    errU=getL2Error(uu,ureal,delta,m);
    errU=sqrt(errU);
    
    
   
    
    
    
     
    
    figure 
    for i=1:n
        sx=[(i-1)*h,i*h];
        sy=[pp(2*i-1),pp(2*i)];
        o1=line(sx,sy,'LineWidth',2);
        hold on
    end
    hold on
    o2=plot(0:h:1,preal(0:h:1),'Color',[1 0.5 0],'LineWidth',1);
    legend([o1,o2],{'numerical-p', 'real-p'})
    xlabel('x')
    ylabel('function value')
    title('NumSol-P')
    
    figure 
    for i=1:n
        sx=[(i-1)*h,i*h];
        sy=[uu(i),uu(i)];
        o1=line(sx,sy,'LineWidth',2);
        hold on
    end
    hold on
    o2=plot(0:h:1,ureal(0:h:1),'Color',[1 0.5 0],'LineWidth',1);
    legend([o1,o2],{'numerical-u', 'real-u'})
    xlabel('x')
    ylabel('function value')
    title('NumSol-U')
    
    
    