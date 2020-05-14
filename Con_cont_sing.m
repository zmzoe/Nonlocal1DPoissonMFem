%Main
% 
%    p-Gu = g
%    -G'p = f-Gg
%
% continuous P1- P0
%
%  u=cos(2 pi x) 
%remove the spcious mode (add \int_Omega uv)
%%-------------------------%%
%% Parameter

delta  =   0.25;
m      =   4;
h      =   delta/m;
n      =   1/h;
k      =   1;
beta   =   1;

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

DD = zeros(n-1,n-1);

for i=1:n-1
    for j = 1:n-1
         DD(i,j) =0.0001*jump(Alpha_p0(j,:), delta, m, i )^2;       
    end
end

MM=zeros(2*n-2,2*n-2);
MM(1:n-1,1:n-1)=A;
MM(1:n-1,n:2*n-2)=-B;
MM(n:2*n-2,1:n-1)=-B';
MM(n:2*n-2,n:2*n-2)=DD;



    
     ureal = @(x) cos(2*pi*k*x);  
     preal = @(x) -2*k*pi*sin(2*pi*k*x);
    
    
     f     = @(x) (2*cos(2*pi*k*x)*(cos(2*pi*delta*k) - 2*delta^2*k^2*pi^2 + 2*delta*k*pi*sin(2*pi*delta*k) - 1))/(delta^4*k^2*pi^2);
 
    
    
   
     g     = @(x)(2*(delta*cos(2*pi*k*x) - (sin(2*pi*k*x) + sin(2*k*pi*(delta - x)))/(2*k*pi)))/delta^2 - 2*k*pi*sin(2*pi*k*x);
    
    
     HH     = @(x) -(2*(cos(2*pi*k*(delta + x)) - cos(2*pi*k*x) + 2*delta*k*pi*sin(2*pi*k*x)))/delta^2;
    
    gq  =  zeros(n-1,1);
    for j=1:n-1
        gq(j,1)=Gauss_anybase(Alpha(j,:),g,h); %(g,q)
    end
    
    HHv  =  zeros(n-1,1);
    for j=1:n-1
        HHv(j,1)=Gauss_anybase(Alpha_p0(j,:),HH,h); %(HH,v)
    end
    
    F     =   [gq;HHv]; 
    
    
    
    
    
      xx=MM\F;
     %[xx,flag,relres] = gmres(MM,F);  
    
    pp1=xx(2:n-1)-xx(1:n-2);
    ppn=[xx(1);pp1;-xx(n-1);xx(1)];
    %n+1
    
    uu0=xx(n:2*n-2);
    uu1=uu0(2:n-1)-uu0(1:n-2);
    uu=[uu0(1);uu1;-uu0(n-1)];
    %n
   
    
    %%
     errP=getDeltaError_anybase(Alpha,xx,HH,delta,m);
     errP=sqrt(errP);
    
      err_P_L2=getL2Error_local_p(ppn,preal,h);
      err_P_L2=sqrt(err_P_L2);
    
       err_U=getL2Errorlocal(uu,ureal,h);
       err_U=sqrt(err_U);
    
    
    figure
    plot(ppn,'y','LineWidth',2);
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
    
    
    
    
   
    
    
    
    
    
    