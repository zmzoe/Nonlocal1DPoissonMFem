%  Main
% 
%    without compability, which means the bases are phi1, phi2.....
%
%    local (p,q)+(dq,u)=0
%               -(dp,v)=(f,v)
%

% discontinuous p1 p0
% Good slop but there is a fix distance between the real solution for p 
%% Parameter

n     =   10;
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


for i = 1:2*n
    for j = 1:n
         B(i,j) = Ponelocal(Alpha(i,:),Alpha_p0(j,:),n );       
    end
end






    ureal = @(x) sin(2*pi*x);  
    preal = @(x) 2*pi*cos(2*pi*x);   
    
    f     = @(x) 4*pi^2*sin(2*pi*x);
   
    

    F     =   Fgauss(h,f); 
      
    
    
        
    FF=[zeros(2*n,1);F];
    
    
    
    MM=zeros(3*n,3*n);
 
    MM(1:2*n,1:2*n)=A(1:2*n,1:2*n);
    MM(2*n+1:3*n,1:2*n)=-B(1:2*n,1:n)';
    MM(1:2*n,2*n+1:3*n)=B(1:2*n,1:n);
    
    
    
    
    
    
    
    % xx=MM\FF;
    [xx,flag,relres] = gmres(MM,FF);
     c1 = sum(xx(1:2*n))/(2*n);
     pp = xx(1:2*n)-c1*ones(2*n,1);

     uu=xx(2*n+1:3*n);
   
    figure   
    plot(pp(1:2:2*n-1));hold on
    plot(preal(0:h:1-h));
    legend('pnum','preal')
    
    
    figure 
   for i=1:n
        sx=[(i-1)*h,i*h];
        sy=[pp(2*i-1),pp(2*i)];
        sz=[preal(sx(1)),preal(sx(2))];
        sw=sz-sy;
        o1=line(sx,sy,'LineWidth',2);
        hold on
        o2=line(sx,sz,'Color',[1 0.5 0],'LineWidth',2);
        hold on
        o3=line(sx,sw,'Color',[1 0.8 1],'LineWidth',2);
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
        sz=[ureal(sx(1)),ureal(sx(2))];
        o1=line(sx,sy,'LineWidth',2);
        hold on
        o2=line(sx,sz,'Color',[1 0.5 0],'LineWidth',2);
        hold on
        o3=line(sx,sz-sy,'Color',[1 0.8 1],'LineWidth',2);
        hold on
        o4=plot((sx(1)+sx(2))/2,(sy(1)+sy(2))/2,'k*');
        hold on
    end
    legend([o1,o2,o3,o4],{'numericalU', 'realU','errU','midpoint'})
    xlabel('x')
    ylabel('function value')
    title('NumSolU')