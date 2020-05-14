%Main
% 

%% Parameter

delta  =        0.5;
m      =         10;
h      =   delta/m;
n      =   1/h;


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

B = -B; % Variontional formula
    
    
    MM = zeros(3*n,3*n);
 
    MM(1:2*n,1:2*n) = A;
    MM(2*n+1:3*n,1:2*n) = B';
    MM(1:2*n,2*n+1:3*n) = B;

%% Real solution and RHS

    ureal = @(x) sin(2*pi*x);  
   preal = @(x) -(2*((cos(2*pi*x) - cos(2*pi*(delta - x)))/(2*pi) + delta*sin(2*pi*x)))/delta^2;    
    f     = @(x) (4*sin(2*pi*x)*(cos(2*pi*delta)/2 - delta^2*pi^2 + delta*pi*sin(2*pi*delta) - 1/2))/(delta^4*pi^2);
   
   
    F  =   Fgauss(h,f);  
    c3 =   sum(F)/n;
    F  =   F-c3*ones(n,1);
    
    %Fhere = F(1:n-1,1)-F(2:n,1);
    FF    = [zeros(2*n,1);F];
    
    
    
   %xx =  MM\FF;
    %[xx,flag,relres] = lsqr(MM,FF);
    [xx,flag,relres] = gmres(MM,FF);
%% Error estimate
     c1 = sum(xx(1:2*n))/(2*n);
     pp = xx(1:2*n)-c1*ones(2*n,1);
     
     pp = [pp;pp(1)];
     
     c2 = sum(xx(2*n+1:3*n))/n;
     uu = xx(2*n+1:3*n)-c2*ones(n,1);
     
   figure
   plot(xx);
    
    %{
    figure
    plot(h/2:h:1-h/2,uu,'y','LineWidth',2);
    hold on
    plot(h/2:h:1-h/2,ureal(h/2:h:1-h/2),'k');
    legend('numerical-u', 'real-u')
    xlabel('x')
    ylabel('function value')
    title('NumSol-U')
    %}
    
    figure
    plot(0:h:1,pp(1:2:2*n+1),'y','LineWidth',2);
    hold on
    plot(0:h:1,preal(0:h:1),'k');
    legend('numerical-p', 'real-p')
    xlabel('x')
    ylabel('function value')
    title('NumSol-P')
    

    
    
  %% err
    
    errP=getDeltaError_dg(Alpha,xx,f,delta,m);
    errP=sqrt(errP);
    %{
    errU=getL2Error(uu,ureal,delta,m);
    errU=sqrt(errU);
   
   %}
    

