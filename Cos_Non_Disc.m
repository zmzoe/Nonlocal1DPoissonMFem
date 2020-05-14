%   Main
% 
%    p-Gu = g
%    -G'p = f-Gg
%
% 
% ureal=cos(2 pi x)
%% Parameter

delta  =   0.5;
m      =   20;
h      =   delta/m;
n      =   1/h;
k      =   1;

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
         A(i,j) = innerproduct(Alpha(i,:),Alpha(j,:), delta, m );       
    end
end


B = zeros(2*n-1,n-1);

for i = 1:2*n-1
    for j = 1:n-1
         B(i,j) = Pone(Alpha(i,:),Alpha_p0(j,:), delta, m );       
    end
end



MM=zeros(3*n-2,3*n-2);
MM(1:2*n-1,1:2*n-1)=A;
MM(1:2*n-1,2*n:3*n-2)=-B;
MM(2*n:3*n-2,1:2*n-1)=-B';






   %ureal = @(x) cos(2*pi*x);  
    %preal = @(x) -2*pi*sin(2*pi*x);   
     ureal = @(x) cos(2*pi*k*x);  
     preal = @(x) -2*k*pi*sin(2*pi*k*x);
    
    %f     = @(x) (2*cos(2*pi*x)*(cos(2*pi*delta) - 2*delta^2*pi^2 + 2*delta*pi*sin(2*pi*delta) - 1))/(delta^4*pi^2);
     f     = @(x) (2*cos(2*pi*k*x)*(cos(2*pi*delta*k) - 2*delta^2*k^2*pi^2 + 2*delta*k*pi*sin(2*pi*delta*k) - 1))/(delta^4*k^2*pi^2);
 
    
    
    %g     = @(x) (2*(delta*cos(2*pi*x) - (sin(2*pi*x) + sin(2*pi*(delta - x)))/(2*pi)))/delta^2 - 2*pi*sin(2*pi*x);  
     g     = @(x)(2*(delta*cos(2*pi*k*x) - (sin(2*pi*k*x) + sin(2*k*pi*(delta - x)))/(2*k*pi)))/delta^2 - 2*k*pi*sin(2*pi*k*x);
    
    %HH     = @(x) -(2*(cos(2*pi*(delta + x)) - cos(2*pi*x) + 2*delta*pi*sin(2*pi*x)))/delta^2;
     HH     = @(x) -(2*(cos(2*pi*k*(delta + x)) - cos(2*pi*k*x) + 2*delta*k*pi*sin(2*pi*k*x)))/delta^2;

    gq  =  zeros(2*n-1,1);
    for j=1:2*n-1
        gq(j,1)=Gauss_anybase(Alpha(j,:),g,h); %(g,q)
    end
    
    HHv  =  zeros(n-1,1);
    for j=1:n-1
        HHv(j,1)=Gauss_anybase(Alpha_p0(j,:),HH,h); %(HH,v)
    end
    
    F     =   [gq;HHv]; 



    
    
    xx=MM\F;
       
    
    
    pp1=xx(2:2*n-1)-xx(1:2*n-2);
    pp=[xx(1);pp1;-xx(2*n-1)]; %pp: 2*n
   
    ppn=[pp;pp(1)];
   
    uu0=xx(2*n:3*n-2);
    uu1=uu0(2:n-1)-uu0(1:n-2);
    uu=[uu0(1);uu1;-uu0(n-1)]; %uu: n
   
    
    
     errP=getDeltaError_anybase_dg(Alpha,xx,HH,delta,m);
     errP=sqrt(errP);
    
      err_P_L2=getL2Error_local_p_dg(ppn,preal,h);
      err_P_L2=sqrt(err_P_L2);
    
    
    
    errU=getL2Error(uu,ureal,delta,m);
    errU=sqrt(errU);
    
    
   
   %% Figure 1------- P ---------P_h---------- P-P_h
    
    
   kkangle=zeros(n,1);  
    
    figure 
   for i=1:n
        sx=[(i-1)*h,i*h];
        sy=[pp(2*i-1),pp(2*i)];
        sz=[preal(sx(1)),preal(sx(2))];
        sw=sz-sy;
        kkangle(i)=abs((sz(2)-sz(1))-(sy(2)-sy(1)))/abs(h^2+(sz(2)-sz(1))*(sy(2)-sy(1)));
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
    
    
    %% Figure 2 =---------p_ph_even-------p_ph_odd------angle------
    errp_ph_even  =  preal(h:h:1)-pp(2:2:2*n)';
    errp_ph_odd   =  preal(0:h:1-h)-pp(1:2:2*n-1)';
    
    figure
    plot(0:h:1-h,errp_ph_odd);
    hold on
    plot(h:h:1,errp_ph_even,'Color',[0.5 0.5 0.5]);
    hold on
    plot(h/2:h:1-h/2,kkangle,'Color',[1 0.8 1],'LineWidth',1.5);
    grid on
    
    legend('pphodd','ppheven','angle')
    xlabel('x')
    ylabel('function value')
    title('VariousErrofP')
   
    
   %% Figure 3 =---------u_uh-----------------
      
    
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
    
    
    
    
   
    
    
