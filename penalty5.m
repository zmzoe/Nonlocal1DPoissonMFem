% 
%  极限到local情况 加稳定项 impluse c0




%     (p,q)+(dq/dx,u)+ sum_i \lambda_i [q]_i=0
%     (dp/dx, v)=(f,v)-gamma h^rho
%     sum_i \nv_i [p]_i=gamma h^rho
% 
% p and u are both simple trigonometric functions
% discontinuous p1 p0

%% Parameter

n     =   20;
h     =   1/n;
gamma =  0.1;
rho   =   2;


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
         A(i,j) = innerproductlocal(Alpha(i,:),Alpha(j,:), n);       
    end
end


B = zeros(2*n,n);


for i = 1:2*n
    for j = 1:n
         B(i,j) = Ponelocal(Alpha(i,:),Alpha_p0(j,:),n );       
    end
end






    ureal = @(x) sin(2*pi*x);  
    preal = @(x) 2*pi*cos(2*pi*x);   
    
    f     = @(x) -4*pi^2*sin(2*pi*x); % laplace u=f
   
    

    F     =   Fgauss(h,f); 
    F     =   F-gamma*h^rho*ones(n,1) ;
      
    
    
        
   
    
    
    
    MM=zeros(3*n,3*n);
 
    MM(1:2*n,1:2*n)=A(1:2*n,1:2*n);
    MM(2*n+1:3*n,1:2*n)=B(1:2*n,1:n)';
    MM(1:2*n,2*n+1:3*n)=B(1:2*n,1:n);
    
    
    CC=zeros(2*n,n);
    for i = 1:n-1
        CC(2*i-1,i) = -1;
        CC(2*i,i+1) =  1;
    end
    
    CC(2*n-1,n) = -1;
    CC(2*n,1)   =  1;
    
    NN=zeros(4*n,4*n);
    NN(1:3*n,1:3*n)=MM;
    NN(1:2*n,3*n+1:4*n)=CC;
    NN(3*n+1:4*n,1:2*n)=CC';
    
    
    FF=[zeros(2*n,1);F;gamma*h^rho*ones(n,1)];
    
    %xx=NN\FF;
    [xx,flag,relres] = gmres(NN,FF);
     c1 = sum(xx(1:2*n))/(2*n);
     pp = xx(1:2*n)-c1*ones(2*n,1);

     uu=xx(2*n+1:3*n);
 
     
     
     
     

    
    
    figure 
   for i=1:n
        sx=[(i-1)*h,i*h];
        sy=[pp(2*i-1),pp(2*i)];
        sz=[preal(sx(1)),preal(sx(2))];
        sw=sz-sy;
        o2=line(sx,sz,'Color',[1,1,0],'LineWidth',4);
        hold on
        o1=line(sx,sy,'Color',[0 0 0],'LineWidth',1);
        hold on
        o3=line(sx,sw,'Color',[1 0.8 1],'LineWidth',1);
        hold on
        o4=plot((sx(1)+sx(2))/2,(sy(1)+sy(2))/2,'k*');
        hold on
   end     
    legend([o1,o2,o3,o4],{'numericalP', 'realP','err','numPmidpoint'})
    xlabel('x')
    ylabel('function value')
    title('NumSolP')
    
    figure 
    for i=1:n
        sx=[(i-1)*h,i*h];
        sy=[uu(i),uu(i)];
        sw=[ureal(sx(1)),ureal(sx(2))]-sy;      
        o1=line(sx,sy,'Color',[0,0,0],'LineWidth',1);
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
   
    figure 
    plot(-xx(3*n+1:4*n))
    hold on
    plot(uu)
    legend('-lambda','u')