% Main
% real u = sin (2*pi*x)
% real p = Gu 
%      p =-(2*((cos(2*pi*x) - cos(2*pi*(delta - x)))/(2*pi) + delta*sin(2*pi*x)))/delta^2
%      
% RHS  f = -G^*G u
%      f =(4*sin(2*pi*x)*(cos(2*pi*delta)/2 - delta^2*pi^2 + delta*pi*sin(2*pi*delta) - 1/2))/(delta^4*pi^2)
%      
%
%% -------------------------------------------------------------------------

    delta = 0.2;
   
    m = 5;
    
    h = delta/m;
    n = 1/h;
    
    %% A ------------------------------------------------------------------
    A = zeros(n,n);
    A(1,1) = 2*h/3;
    A(1,n) = h/6;
    A(n,1) = h/6;
    
   
    
    
    for i = 2:n
        A(i,i) = 2*h/3;
    end
    
    for i = 1:n-1
        A(i,i+1) = h/6;   
        A(i+1,i) = h/6;
    end
    
    
    B = sparse(n,n);
    
    vv = zeros(m+2,1);
    for i = 1:n
    I = i*ones(1,m+2);
   
    
    if i<=m+1
        J = [n-m+i-1:n,1:i];
    else
        J = i-m-1:i;
    end
    
    
    if m == 1
        vv(1) = h^2/6;
        vv(2) = h^2/6;
        vv(3) = -h^2/3;
    else
        vv(1) = h^2/6;
        vv(2) = h^2*5/6;
    for j=3:m
        vv(j) = h^2;
    end
    
    vv(m+1) = 5*h^2/6-m*h^2/2;
    vv(m+2) = h^2/6-m*h^2/2;
   
    end
    
    B = B+sparse(I,J,vv,n,n);
   
    
    end  
      
    
    B = -2/(delta^2)*B;  
    B = full(B);
    
  % Real solution
    ureal = @(x) sin(2*pi*x);  
    preal = @(x) -(2*((cos(2*pi*x) - cos(2*pi*(delta - x)))/(2*pi) + delta*sin(2*pi*x)))/delta^2;    
    f=@(x) (4*sin(2*pi*x)*(cos(2*pi*delta)/2 - delta^2*pi^2 + delta*pi*sin(2*pi*delta) - 1/2))/(delta^4*pi^2);
           
    F=Fgauss(h,f);          
    FF = [zeros(n,1);F];
    
    
    MM = zeros(2*n,2*n);
 
    MM(1:n,1:n) = A;
    MM(n+1:2*n,1:n) = B';
    MM(1:n,n+1:2*n) = B;
    
    
    %xx=MM\FF;
    [xx,flag,relres] = minres(MM,FF); 
    
    pp=[xx(1:n);xx(1)];
    uu=xx(n+1:2*n);
        
    figure
    plot(xx)
    
    
    figure
    plot(0:h:1,pp,'y','LineWidth',2);
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

   
    %% err
    
    errP=getDeltaError(xx,f,delta,m);
    errP=sqrt(errP);
   
    errU=getL2Error(uu,ureal,delta,m);
    errU=sqrt(errU);
