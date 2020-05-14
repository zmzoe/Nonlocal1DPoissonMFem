%   Main
% 
%         [(1+h^2)(p,q)+mu_2(dp/dx,dq/dx)]+(1+h^2)(dq/dx,u)=-mu_2(f,dq/dx)
%         for any q
%         [(dp/dx,v)+mu_1 h^2(dp/dx,dv/dx)]-mu_1 h^2(du/dx,dv/dx)=-(f,v)
%         for any v
%        
% . when mu_2=mu_1=0, the variational formular is the same as the 1D local
% case
% 
% 
% discontinuous p1 p0
% base functions are phi1.....phi2n 

%% Parameter

n     =   20;
h     =   1/n;
t     =   -0.5;
mu    =   1;
mu1   =   0.5;
mu2   =   1;


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
         A(i,j) = (1+mu*h^2*t)*innerproductlocal(Alpha(i,:),Alpha(j,:), n)+(mu2*mu*h^2)*innerDPproductlocal(Alpha(i,:),Alpha(j,:), n);       
    end
end



B = zeros(2*n,n);

for i = 1:2*n
    for j = 1:n
         B(i,j) = Ponelocal(Alpha(i,:),Alpha_p0(j,:),n );       
    end
end


B  =  (1+mu*h^2*t)*B;


%%


    ureal = @(x) sin(2*pi*x);  
    preal = @(x) 2*pi*cos(2*pi*x);   
    
    f     = @(x) 4*pi^2*sin(2*pi*x);
   
    

    F     =   Fgauss(h,f); 
      
    
    FF1 = zeros(2*n,1);
    for k = 1:n
        FF1(2*k-1,1) =  1/h*F(k);
        FF1(2*k,1)   = -1/h*F(k);
    end
    
    FF1 = mu2*mu*h^2*FF1;
        
    FF=[FF1;F];
    
    
    C = zeros(2*n,n);

for i = 1:2*n
    for j = 1:n
         C(i,j) = Ponelocal(Alpha(i,:),Alpha_p0(j,:),n );       
    end
end

    
    
    
    MM=zeros(3*n,3*n);
 
    MM(1:2*n,1:2*n)=A(1:2*n,1:2*n);
    MM(2*n+1:3*n,1:2*n)=-C(1:2*n,1:n)';
    MM(1:2*n,2*n+1:3*n)=B(1:2*n,1:n);
    
    
    
    
    
    
    
     %xx=MM\FF;
    [xx,flag,relres] = gmres(MM,FF);
     c1 = sum(xx(1:2*n))/(2*n);
     pp = xx(1:2*n)-c1*ones(2*n,1);

     uu=xx(2*n+1:3*n);
   
     figure
     plot(uu)
     
     
   
    
    
    
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