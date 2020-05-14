%%-------------------------%%
%% Parameter

delta  =   0.5;
m      =   5;
h      =   delta/m;
n      =   1/h;

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

dd=sum(A,2);

Amasslump=diag(dd);
B = zeros(n-1,n-1);

for i = 1:n-1
    for j = 1:n-1
         B(i,j) = Pone(Alpha(i,:),Alpha_p0(j,:), delta, m );       
    end
end



MM=zeros(2*n-2,2*n-2);
MM(1:n-1,1:n-1)=A;
MM(1:n-1,n:2*n-2)=-B;
MM(n:2*n-2,1:n-1)=-B';

%%
ureal={@(x) (-x+1/4);
       @(x) (x-3/4)};

   
% (0,delta)
% (delta,1/2)
% (1/2,1/2+delta)
% (1/2+delta,1)
preal={ @(x) (2*(delta.*(x - 1/4) - (x.*(2.*x - 1))/4 + ((delta - x).*(2.*x - 2.*delta + 1))/4))/delta^2;
        @(x) 1;
        @(x) (2.*(((2.*x - 1).*(x - 1))/4 - delta.*(x - 3/4) + ((delta - x).*(2.*delta - 2.*x + 1))/4))/delta^2;
        @(x) -1};

%delta<=1/4 ÐèÒª·ÖÁù¶Î
%(0,delta) 
%(delta,1/2-delta) 
%(1/2-delta,1/2) 
%(1/2,1/2+delta) 
%(1/2+delta,1-delta)

   if delta>1/4 && delta<=1/2  
       f={@(x) (-8/(3*delta^4))*x.^3 + (2/delta^4)*x.^2 + (-(1/delta^2 - 8)/delta^2).*x + (2*delta - (2*(delta/2 + (delta*(2*delta - 1))/2))/delta + 1/(6*delta^2) - 2)/delta^2;
          @(x) (-4/(3*delta^4))*x.^3 + (2/delta^4)*x.^2 + (-(1/delta^2 - 4)/delta^2).*x + ((8*delta)/3 + 1/(6*delta^2) - 2)/delta^2;
          @(x) (-4/(3*delta^4))*x.^3 + (2/delta^4)*x.^2 + (-(1/delta^2 - 4)/delta^2).*x + ((8*delta)/3 + 1/(6*delta^2) - 2)/delta^2;
          @(x) (8/(3*delta^4))*x.^3 + (-((2*(2*delta + 1))/delta^2 - 4/delta + 4/delta^2)/delta^2)*x.^2 + (((2*(2*delta + 1/2))/delta^2 - (2*(4*delta + 2))/delta + 4/delta^2)/delta^2)*x - (2*delta - (2*((3*delta)/2 + (delta*(2*delta + 1))/2 + 1/2))/delta + (2*(delta/2 + 1/12))/delta^2 + 4/(3*delta^2) - 2)/delta^2;
          @(x) (4/(3*delta^4))*x.^3 + (-4/delta^4)*x.^2 + ((4/delta^2 - 4)/delta^2).*x - ((8*delta)/3 + 4/(3*delta^2) - 4)/delta^2;
          @(x) (4/(3*delta^4))*x.^3 + (-4/delta^4)*x.^2 + ((4/delta^2 - 4)/delta^2).*x - ((8*delta)/3 + 4/(3*delta^2) - 4)/delta^2};

   elseif delta>0 && delta<=1/4
           f={@(x) (-4/(3*delta^4))*x.^3 + (4/delta^2).*x - ((2*delta)/3 + (2*(delta/2 + (delta*(2*delta - 1))/2))/delta)/delta^2;
              @(x) 0;
              @(x) (-4/(3*delta^4))*x.^3 + (2/delta^4)*x.^2 + (-(1/delta^2 - 4)/delta^2).*x + ((8*delta)/3 + 1/(6*delta^2) - 2)/delta^2; 
              @(x) (4/(3*delta^4))*x.^3 + (-((2*(2*delta + 1))/delta^2 - 4/delta)/delta^2).*x.^2 + (((2*(2*delta + 1/2))/delta^2 - (2*(4*delta + 2))/delta + 4)/delta^2)*x + ((2*delta)/3 + (2*((3*delta)/2 + (delta*(2*delta + 1))/2 + 1/2))/delta - (2*(delta/2 + 1/12))/delta^2 - 2)/delta^2;
              @(x) 0;
              @(x) (4/(3*delta^4)).*x.^3 + (-4/delta^4).*x.^2 + ((4/delta^2 - 4)/delta^2).*x - ((8*delta)/3 + 4/(3*delta^2) - 4)/delta^2};
 
     
     
   end
               
     
%%   
     
    
    HHv  =  zeros(n-1,1);
    for j=1:n-1
        HHv(j,1)=Gauss_anybasePoly(Alpha_p0(j,:),f,h,delta); %(HH,v)
    end
    
    F     =   [zeros(n-1,1);HHv]; 
    
    
    
    
    
      xx=MM\F;
      
     %[xx,flag,relres] = gmres(MM,F);  
    
    pp1=xx(2:n-1)-xx(1:n-2);
    ppn=[xx(1);pp1;-xx(n-1);xx(1)];
    %n+1
    
    uu0=xx(n:2*n-2);
    uu1=uu0(2:n-1)-uu0(1:n-2);
    uu=[uu0(1);uu1;-uu0(n-1)];
    %n
   
    
    %%error
    
     errP=getDeltaError_anybasePoly(Alpha,xx,f,delta,m);
     errP=sqrt(errP);
     a1=errP;
    
     err_P_L2=getL2Error_local_pPOLY(ppn,preal,h,delta);
     err_P_L2=sqrt(err_P_L2);
     a2=err_P_L2;
    
     
     err_U=getL2ErrorlocalPOLY(uu,ureal,h);
     err_U=sqrt(err_U);
     a3=err_U;
     
    
   %{
    
   figure 
    for i=1:n
        sx=[(i-1)*h,i*h];
        sy=[uu(i),uu(i)];
        o1=line(sx,sy,'LineWidth',2);      
        hold on
    end
    hold on
    o2=plot(0:h:1/2,ureal{1}(0:h:1/2),'Color',[1 0.5 0],'LineWidth',1);
    hold on
    o3=plot(1/2:h:1,ureal{2}(1/2:h:1),'Color',[1 0.5 0],'LineWidth',1);
    
    legend([o1,o2,o3],{'numerical-u', 'real-u','real-u'})
    
    xlabel('x')
    ylabel('function value')
    title('NumSol-U')
    
    [~,zm]=size(delta:h:1/2);
    ppreal=zeros(zm,2);
    for i=1:zm
        ppreal(i,1)=preal{2}(delta+(i-1)*h);
        ppreal(i,2)=preal{4}(delta+(i-1)*h);
    end
   
    figure
    plot(0:h:1,ppn,'k','LineWidth',2);
    hold on
    plot(0:h:delta,preal{1}(0:h:delta),'r','LineWidth',1);
    hold on
    plot(delta:h:1/2,ppreal(:,1),'r','LineWidth',1);
    hold on
    plot(1/2:h:1/2+delta,preal{3}(1/2:h:1/2+delta),'r','LineWidth',1);
    hold on
    plot(1/2+delta:h:1,ppreal(:,2),'r','LineWidth',1);
    
    legend({'numerical-p', 'real-p','real-p'})
    xlabel('x')
    ylabel('function value')
    title('NumSol-P')
   
    
    
    figure
    plot(0:h:delta,f{1}(0:h:delta),'r','LineWidth',2);
    hold on
    plot(delta:h:1/2-delta,f{2}(delta:h:1/2-delta),'r','LineWidth',2);
    hold on
    plot(1/2-delta:h:1/2,f{3}(1/2-delta:h:1/2),'r','LineWidth',2);
    hold on
    plot(1/2:h:1/2+delta,f{4}(1/2:h:1/2+delta),'r','LineWidth',2);
    hold on
    plot(1/2+delta:h:1-delta,f{5}(1/2+delta:h:1-delta),'r','LineWidth',2);
    hold on
    plot(1-delta:h:1,f{6}(1-delta:h:1),'r','LineWidth',2);
    
    xlabel('x')
    ylabel('function value')
    title('real f')
    
    figure
    plot(0:h:delta,preal{1}(0:h:delta),'r','LineWidth',1);
    hold on
    plot(delta:h:1/2,ppreal(:,1),'r','LineWidth',1);
    hold on
    plot(1/2:h:1/2+delta,preal{3}(1/2:h:1/2+delta),'r','LineWidth',1);
    hold on
    plot(1/2+delta:h:1,ppreal(:,2),'r','LineWidth',1);
    
    legend({'numerical-p', 'real-p','real-p'})
    xlabel('x')
    ylabel('function value')
    title('NumSol-P')
   
    figure
    o2=plot(0:h:1/2,ureal{1}(0:h:1/2),'Color',[1 0.5 0],'LineWidth',1);
    hold on
    o3=plot(1/2:h:1,ureal{2}(1/2:h:1),'Color',[1 0.5 0],'LineWidth',1);
    
   %} 