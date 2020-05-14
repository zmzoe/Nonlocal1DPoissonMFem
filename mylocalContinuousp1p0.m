%   Main
% 
%   
%
% 
% continuous p1 p0

%% Parameter

n     =   4;
h     =   1/n;


%% basis 
%% P1
Alpha = zeros(n, 2*n);

for i = 2:n
    Alpha(i,2*(i-1)) = 1;
    Alpha(i,2*i-1)   = 1;
end

Alpha(1,1)  =  1;
Alpha(1,2*n)=  1;

% P0
Alpha_p0 = zeros(n, 2*n);

for i = 1:n
    Alpha_p0(i,2*i-1)   = 1;
    Alpha_p0(i,2*i)     = 1;
 
end


%% Linear system

A = zeros(n,n);

for i = 1:n
    for j = 1:n
         A(i,j) = innerproductlocal(Alpha(i,:),Alpha(j,:), n);       
    end
end

Astiffness = zeros(n,n);

for i = 1:n
    for j = 1:n
         Astiffness(i,j) = Stifnesslocal(Alpha(i,:),Alpha(j,:), n);   
    end
end


B = zeros(n,n);


for i = 1:n
    for j = 1:n
         B(i,j) = Ponelocal(Alpha(i,:),Alpha_p0(j,:),n );       
    end
end






    ureal = @(x) sin(2*pi*x);  
    preal = @(x) 2*pi*cos(2*pi*x);   
    
    f     = @(x) 4*pi^2*sin(2*pi*x);
   
    

    F     =   Fgauss(h,f); 
    
    FF    =   [zeros(n,1);F];
    
    
    MM=zeros(2*n,2*n);
 
    MM(1:n,1:n)=A;
    MM(n+1:2*n,1:n)=-B';
    MM(1:n,n+1:2*n)=B;
    
    
    
    [xx,flag,relres] = gmres(MM,FF);
    
    pp=xx(1:n);
    
    pp1=[pp;xx(1)];
    
    uu=xx(n+1:2*n);
    
    
    err_P_H1=getH1Error_local_p(pp1,f,h);
    err_P_H1=sqrt(err_P_H1);
    
    err_P_L2=getL2Error_local_p(pp1,preal,h);
    err_P_L2=sqrt(err_P_L2);
    
    err_U=getL2Errorlocal(uu,ureal,h);
    err_U=sqrt(err_U);
    
    
    p_real=preal(linspace(0,1-h,n)');
    
    
    coe=p_real-pp;

    
err=sqrt(coe'*A*coe);

    
    figure
    plot(pp1,'y','LineWidth',2);
    hold on
    plot(preal(0:h:1));
    legend('numP','realP')
    
    
   figure 
    for i=1:n
        sx=[(i-1)*h,i*h];
        sy=[uu(i),uu(i)];
        sw=[ureal(sx(1)),ureal(sx(2))]-sy;
        o1=line(sx,sy,'LineWidth',2);
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
    
    
    
    
   
    
    
    