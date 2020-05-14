

%% Mass lumpping
%------------------------------------------
% 
%    p-Gu = 0
%    -G'p = f
%------------------------------------------
% continuous P1- P0
%------------------------------------------
%  u=cos(2 pi x) 
%------------------------------------------
% base function consists of hat functions
%%-----------------------------------------


%% Parameter

delta  =   0.5;
m      =   5;
h      =   delta/m;
n      =   1/h;
k      =   2;

result = cell(4,1);


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


MM=zeros(2*n,2*n);
MM(1:n,1:n)=Amasslump;
MM(1:n,n+1:2*n)=-B;
MM(n+1:2*n,1:n)=-B';

%{
%% Real-solution
ureal = @(x) cos(2*pi*k*x);  
preal = @(x) -(2*(delta*cos(2*pi*k*x) - (sin(2*pi*k*x) + sin(2*k*pi*(delta - x)))/(2*k*pi)))/delta^2;

f = @(x) (2*cos(2*pi*k*x)*(cos(2*pi*delta*k) - 2*delta^2*k^2*pi^2 + 2*delta*k*pi*sin(2*pi*delta*k) - 1))/(delta^4*k^2*pi^2);
%}

ureal={@(x) (-x+1/4);
       @(x) (x-3/4)};
   
preal={ @(x) (2*(delta.*(x - 1/4) - (x.*(2.*x - 1))/4 + ((delta - x).*(2.*x - 2.*delta + 1))/4))/delta^2;
        @(x) 1;
        @(x) (2.*(((2.*x - 1).*(x - 1))/4 - delta.*(x - 3/4) + ((delta - x).*(2.*delta - 2.*x + 1))/4))/delta^2;
        @(x) -1};
    
    
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
%{
 HHv  =  zeros(n,1);
    for j=1:n
        HHv(j,1)=Gauss_anybase(Alpha_p0(j,:),f,h); %(f,v)
    end
    %}
   
   
    HHv  =  zeros(n,1);
    for j=1:n
        HHv(j,1)=Gauss_anybasePoly(Alpha_p0(j,:),f,h,delta); %(HH,v)
    end
    
    
%% Matrix system 

F     =   [zeros(n,1);HHv];

%xx=MM\F;

[xx,flag,relres] = gmres(MM,F);

ppn=[xx(1:n);xx(1)];
uu=xx(n+1:2*n);



    
     err_P_L2=getL2Error_local_pPOLY(ppn,preal,h,delta);
     err_P_L2=sqrt(err_P_L2);
     a1=err_P_L2;
     
     err_U=getL2ErrorlocalPOLY(uu,ureal,h);
     err_U=sqrt(err_U);
     a2=err_U;

%{

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
    %}

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
    