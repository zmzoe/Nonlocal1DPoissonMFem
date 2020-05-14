% Main
% Rank discontinuous A B

%% Parameter

delta  =   0.5;
m      =   10;
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


     preal = @(x) sin(2*pi*x);  
     f =@(x) (2*(delta*sin(2*pi*x) + (cos(2*pi*(delta + x)) - cos(2*pi*x))/(2*pi)))/delta^2;   
   
    F  =   Fgauss(h,f);  
   
   
    
    Fhere = F(1:n-1,1)-F(2:n,1);
    
    
    FF    = [zeros(2*n-1,1);Fhere];
    
   
    
    
    
    
    xx=MM\FF;
    %[xx,flag,relres] = minres(MM,FF);
    
    pp1=xx(2:2*n-1)-xx(1:2*n-2);
    pp=[xx(1);pp1;-xx(2*n-1)];
   
    figure 
    plot(pp);
    
   
    uu0=xx(2*n:3*n-2);
    uu1=uu0(2:n-1)-uu0(1:n-2);
    uu=[uu0(1);uu1;-uu0(n-1)];
   
    
    
    errP=getDeltaError_dgrank(Alpha,xx,f,delta,m);
    errP=sqrt(errP);
    
    
    %errU=getL2Error(uu,ureal,delta,m);
    %errU=sqrt(errU);
    
    
    figure
    plot(xx);
    
    
    figure
    plot(0:h:1-h,pp(1:2:2*n-1),'y','LineWidth',2);
    hold on
    plot(0:h:1,preal(0:h:1),'k');
    legend('numerical-p', 'real-p')
    xlabel('x')
    ylabel('function value')
    title('NumSol-P')
    
