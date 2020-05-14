% Main
% Continuous
% Rank A and B

%% Parameter

delta  =       0.5;
m      =         5;
h      =   delta/m;
n      =       1/h;


%% basis 1 continuous P1- P0
%% P1

Alpha = zeros(n-1, 2*n);

for i = 2:n-1
    Alpha(i,2*(i-1)) =     1;
    Alpha(i,2*i-1)   =     1;
    Alpha(i,2*i)     =    -1;
    Alpha(i,2*i+1)   =    -1;
end

Alpha(1,1)   =      1;
Alpha(1,2)   =     -1;
Alpha(1,3)   =     -1;
Alpha(1,2*n) =      1;


%% P0

Alpha_p0 = zeros(n-1, 2*n);

for i = 1:n-1
    Alpha_p0(i,2*i-1)     =  1;
    Alpha_p0(i,2*i)       =  1;
    Alpha_p0(i,2*i+1)     = -1;
    Alpha_p0(i,2*i+2)     = -1;
 
end

%% Linear system

A = zeros(n-1,n-1);

for i = 1:n-1
    for j = 1:n-1
         A(i,j) = innerproduct(Alpha(i,:),Alpha(j,:), delta, m );       
    end
end


B = zeros(n-1,n-1);

for i = 1:n-1
    for j = 1:n-1
         B(i,j) = Pone(Alpha(i,:),Alpha_p0(j,:), delta, m );       
    end
end

%aa   =   rank(A);
%aa_1 =   size(A);
%bb   =   rank(B);   
%bb_1 =   size(B);   

B    =   -B;


    MM = zeros(2*n-2,2*n-2);
 
    MM(1:n-1,1:n-1) = A;
    MM(n:2*n-2,1:n-1) = B';
    MM(1:n-1,n:2*n-2) = B;
    

    ureal = @(x) sin(2*pi*x);  
    preal = @(x) -(2*((cos(2*pi*x) - cos(2*pi*(delta - x)))/(2*pi) + delta*sin(2*pi*x)))/delta^2;    
    f     = @(x) (4*sin(2*pi*x)*(cos(2*pi*delta)/2 - delta^2*pi^2 + delta*pi*sin(2*pi*delta) - 1/2))/(delta^4*pi^2);

    F  =   Fgauss(h,f);  
    
    Fhere =  F(1:n-1,1)-F(2:n,1);
    FF    = [zeros(n-1,1);Fhere];  
    
    
    
    
     xx=MM\FF;
   % [xx,flag,relres] = minres(MM,FF);
    
    
    
    
    pp1=xx(2:n-1)-xx(1:n-2);
    ppn=[xx(1);pp1;-xx(n-1)];
    
    figure
    plot(ppn);
    
   
    uu0=xx(n:2*n-2);
    uu1=uu0(2:n-1)-uu0(1:n-2);
    uu=[uu0(1);uu1;-uu0(n-1)];
    figure
    plot(uu);

    
    
    %errP=getDeltaError_dgrank(Alpha,xx,f,delta,m);
    %errP=sqrt(errP);
    
    
    
    %errU=getL2Error(uu,ureal,delta,m);
    %errU=sqrt(errU);
   
    figure
    plot(0:h:1-h,ppn,'y','LineWidth',2);
    hold on
    plot(0:h:1,preal(0:h:1),'k');
    legend('numerical-p', 'real-p')
    xlabel('x')
    ylabel('function value')
    title('NumSol-P')
    
    
    figure
    plot(h/2:h:1-h/2,uu,'y','LineWidth',2);
    hold on
    plot(h/2:h:1-h/2,ureal(h/2:h:1-h/2),'k');
    legend('numerical-u', 'real-u')
    xlabel('x')
    ylabel('function value')
    title('NumSol-U')
    %}