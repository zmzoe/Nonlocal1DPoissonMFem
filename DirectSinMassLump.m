


%% ÷±Ω”º∆À„p=Gu
%  
%  Mass lump 
%%-------------------------%%




%% Parameter

delta  =   0.125/8;
m      =   5;
h      =   delta/m;
n      =   1/h;
k      =   1;


%% basis 1 continuous P1- P0

%% P1

Alpha = zeros(n, 2*n);

for i = 2:n
    Alpha(i,2*(i-1)) =     1;
    Alpha(i,2*i-1)   =     1;
end

Alpha(1,1)   =      1;
Alpha(1,2*n) =      1;


%% P0

Alpha_p0 = zeros(n, 2*n);

for i = 1:n
    Alpha_p0(i,2*i-1)     =  1;
    Alpha_p0(i,2*i)       =  1;
 
end


%% Linear system

A = zeros(n,n);

for i = 1:n
    for j = 1:n
         A(i,j) = innerproduct(Alpha(i,:),Alpha(j,:), delta, m );       
    end
end


dd=sum(A,2);

Amasslump=diag(dd);

B = zeros(n,n);

for i = 1:n
    for j = 1:n
         B(i,j) = Pone(Alpha(i,:),Alpha_p0(j,:), delta, m );       
    end
end
[UU,S,VV] = svd (B);
ww=UU(:,n-1);
wwtest=ones(n,1);
for i=1:2:n-1
    wwtest(i,1)=-1;
end

MM=zeros(2*n,2*n);
MM(1:n,1:n)=A;
MM(1:n,n+1:2*n)=-B;
MM(n+1:2*n,1:n)=-B';

%%

     ureal = @(x) cos(2*pi*k*x).*sin(4*pi*k*x); 
     preal = @(x) -(2*(delta*sin(4*pi*k.*x).*cos(2*pi*k.*x) - (cos(2*k*pi*(delta - x))./4 + cos(6*k*pi.*(delta - x))./12 - cos(2*pi*k.*x)/4 - cos(6*pi*k.*x)/12)./(k*pi)))/delta^2;
    
     f     = @(x)  -(2*(2*sin(4*pi*k*x).*cos(2*pi*k*x) - ((sin(2*pi*k*x).*cos(2*pi*delta*k))/2 - sin(6*pi*k*x)/18 - sin(2*pi*k*x)/2 + (sin(6*pi*k*x).*cos(6*pi*delta*k))/18 + delta*(k*pi*sin(2*pi*k*x).*sin(2*pi*delta*k) + (k*pi*sin(6*pi*k*x).*sin(6*pi*delta*k))/3))/(delta^2*k^2*pi^2)))/delta^2;
     
      
    HHv  =  zeros(n,1);
    for j=1:n
        HHv(j,1)=Gauss_anybase(Alpha_p0(j,:),f,h); %(HH,v)
    end


     
   
    
    F     =   [zeros(n,1);HHv]; 
    
    
    
    
   [xx,flag,relres] = gmres(MM,F);  
    %[xx,flag,relres] = lsqr(MM,F);
    ppn=[xx(1:n);xx(1)];
    uu=xx(n+1:2*n);
    %%
    
    
      err_P_L2=getL2Error_local_p(ppn,preal,h);
      err_P_L2=sqrt(err_P_L2);
      a1=err_P_L2;
      
      
      coe=preal(linspace(0,1-h,n)')-ppn(1:n);
      
      
       err_U=getL2Errorlocal(uu,ureal,h);
       err_U=sqrt(err_U);
       a2=err_U;
  
    figure
    plot(0:h:1,ppn,'LineWidth',2);
    hold on
    plot(0:h:1,preal(0:h:1),'--o','LineWidth',2);
    
    legend('numerical-p', 'real-p')
    xlabel('x')
    ylabel('function value')
    title('NumSol-p')
    
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
    
 
    
    
    
    
    
    