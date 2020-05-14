%   Main
% 
%    with compability, which means the bases are phi1-phi2, phi2-phi3,.....
%
%    local (p,q)+(dq,u)=0
%               -(dp,v)=(f,v)
%
% discontinuous p1 p0

%% Parameter

n     =   20;
h     =   1/n;


%% basis 1 discontinuous P1- P0
%% P1


Alpha = zeros(2*n-1, 2*n);

for i=1:2*n-1
    Alpha(i,i)   =  1;
    Alpha(i,i+1) = -1;
end


%% P0

Alpha_p0 = zeros(n-1, 2*n);

for i = 1:n-1
    Alpha_p0(i,2*i-1)     =  1;
    Alpha_p0(i,2*i)       =  1;
    Alpha_p0(i,2*i+1)     = -1;
    Alpha_p0(i,2*i+2)     = -1;
 
end



%% Linear system

A = zeros(2*n-1,2*n-1);

for i = 1:2*n-1
    for j = 1:2*n-1
         A(i,j) = innerproductlocal(Alpha(i,:),Alpha(j,:), n);       
    end
end


B = zeros(2*n-1,n-1);


for i = 1:2*n-1
    for j = 1:n-1
         B(i,j) = Ponelocal(Alpha(i,:),Alpha_p0(j,:),n );       
    end
end



MM=zeros(3*n-2,3*n-2);
MM(1:2*n-1,1:2*n-1)=A;
MM(1:2*n-1,2*n:3*n-2)=B;
MM(2*n:3*n-2,1:2*n-1)=-B';



    ureal = @(x) sin(2*pi*x);  
    preal = @(x) 2*pi*cos(2*pi*x);   
    
    f     = @(x) 4*pi^2*sin(2*pi*x);
   
    

     F     =   Fgauss(h,f); 
      
     Fhere = F(1:n-1,1)-F(2:n,1);
    
     FF    = [zeros(2*n-1,1);Fhere];
    
        
  
    
    
    
    xx=MM\FF;
    %[xx,flag,relres] = gmres(MM,FF);
    
    pp1=xx(2:2*n-1)-xx(1:2*n-2);
    pp=[xx(1);pp1;-xx(2*n-1)]; %pp: 2*n
   
   
    uu0=xx(2*n:3*n-2);
    uu1=uu0(2:n-1)-uu0(1:n-2);
    uu=[uu0(1);uu1;-uu0(n-1)]; %uu:
   
    figure   
    plot(pp(1:2:2*n-1));
    hold on
    plot(pp(2:2:2*n));
    legend('pnumodd','pnumeven')
    
    figure   
    plot(uu);
    legend('numU')
    
    
    figure 
    for i=1:n
        sx=[(i-1)*h,i*h];
        sy=[pp(2*i-1),pp(2*i)];
        o1=line(sx,sy,'LineWidth',2);
   end     
    
    
    SSWW=zeros(n,2);
    figure 
   for i=1:n
        sx=[(i-1)*h,i*h];
        sy=[pp(2*i-1),pp(2*i)];
        sz=[preal(sx(1)),preal(sx(2))];
        sw=sz-sy;
        SSWW(i,1:2)=sw;
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
    
    NNEE=zeros(n,2);
    figure 
    for i=1:n
        sx=[(i-1)*h,i*h];
        sy=[uu(i),uu(i)];
        sw=[ureal(sx(1)),ureal(sx(2))]-sy;
        NNEE(i,1:2)=[ureal(sx(1))-sy(1),ureal(sx(2))-sy(2)];
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
    
   