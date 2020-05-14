

%---------------------------------------------------%
% Main                                              %
% Test for M9                                       %
% Test for the transformation strategy              %
%---------------------------------------------------%

%% Parameter

delta  =   0.2;
m      =    4;
h      =   delta/m;
n      =   1/h;


%% basis 1 continuous P1- P0
%% P1

Alpha = zeros(n, 2*n);

for i = 2:n
    Alpha(i,2*(i-1)) = 1;
    Alpha(i,2*i-1)   = 1;
end

Alpha(1,1)  =  1;
Alpha(1,2*n)=  1;

%% P0

Alpha_p0 = zeros(n, 2*n);

for i = 1:n
    Alpha_p0(i,2*i-1)   = 1;
    Alpha_p0(i,2*i)     = 1;
 
end

%% Linear system

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

B = -B; % Variontional formula
    
    
    MM = zeros(2*n,2*n);
 
    MM(1:n,1:n) = A;
    MM(n+1:2*n,1:n) = B';
    MM(1:n,n+1:2*n) = B;

%% Real solution and RHS

    ureal = @(x) sin(2*pi*x);  
    preal = @(x) -(2*((cos(2*pi*x) - cos(2*pi*(delta - x)))/(2*pi) + delta*sin(2*pi*x)))/delta^2;    
    f     = @(x) (4*sin(2*pi*x)*(cos(2*pi*delta)/2 - delta^2*pi^2 + delta*pi*sin(2*pi*delta) - 1/2))/(delta^4*pi^2);

    F  =   Fgauss(h,f);  
    
    %Fhere = F(1:n-1,1)-F(2:n,1);
    FF    = [zeros(n,1);F];
    
    
    
   %xx =  MM\FF;
    [xx,flag,relres] = gmres(MM,FF);
%% Error estimate

      uu=xx(n+1:2*n);

   
     errP=getDeltaError_dg(Alpha,xx,f,delta,m);
     errP=sqrt(errP);
    

     errU=getL2Error(uu,ureal,delta,m);
     errU=sqrt(errU);


%% figure
 
figure
plot(xx)
