% 
%  d^2u/dx^2=f

%  ¶ÔÅ¼Íø¸ñ


%     (p,q)+(dq/dx,u)+ sum_i u_i [q]_i=0
%     (dp/dx, v)=(f,v)
%     
% 
% p and u are both simple trigonometric functions
% discontinuous p1 p0

%% Parameter

n     =   20;
h     =   1/n;



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


for k=1:n
    B(2*k-1,k)=-1/2;
    B(2*k,k)=1/2;   
end


for k=1:n-1
    B(2*k,k)=B(2*k,k)-1;
    B(2*k+1,k)=B(2*k+1,k)+1;  
end

B(2*n,n)=B(2*n,n)-1;
B(1,n)=B(1,n)+1;

for k=1:n-1
    B(2*k-1,k+1)=B(2*k-1,k+1)-1/2;
    B(2*k,k+1)=B(2*k,k+1)+1/2;  
end

B(2*n-1,1)=B(2*n-1,1)-1/2;
B(2*n,1)=B(2*n,1)+1/2;

C=zeros(n,2*n);
for i=2:n
    C(i,2*i-3)=-1/2;
    C(i,2*i-1)=-1/2;
    C(i,2*i-2)=1/2;
    C(i,2*i)=1/2;    
end


C(1,1)=-1/2;
C(1,2*n-1)=-1/2;
C(1,2)=1/2;
C(1,2*n)=1/2;

    ureal = @(x) sin(2*pi*x);  
    preal = @(x) 2*pi*cos(2*pi*x);   
    
    f     = @(x) -4*pi^2*sin(2*pi*x); % laplace u=f
   
    

    F     =   FgaussMove(h,f); 
    
      
    
    FF=[zeros(2*n,1);F];
        
   
    
    
    
    MM=zeros(3*n,3*n);
 
    MM(1:2*n,1:2*n)=A(1:2*n,1:2*n);
    MM(2*n+1:3*n,1:2*n)=C(1:n,1:2*n);
    MM(1:2*n,2*n+1:3*n)=B(1:2*n,1:n);
    
    
    
    
   
    
  
   %xx=NN\FF;
    [xx,flag,relres] = gmres(MM,FF);
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
    for i=2:n-1
        sx=[(i-1)*h-h/2,i*h+h/2];
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
    plot(uu)
    
   