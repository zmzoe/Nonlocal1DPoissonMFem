%%  Main
% real u = sin (2*pi*x)
% real p = Gu 
%      p =
%      2/(delta^2)*((cos(2*pi*delta-2*pi*x)-cos(2*pi*x))/pi-delta*sin(2*pi*x))
% RHS  f = -G^*G u
%      f =
%      -2/(delta^2)*(2*sin(2*pi*x)-1/(delta^2*pi^2)*(2*cos(2*pi*delta)*sin(2*pi*x)-2*sin(2*pi*x)+delta*(pi*cos(2*pi*delta)*cos(2*pi*x)-pi*cos(2*pi*x)+3*pi*sin(2*pi*delta)*sin(2*pi*x))));
%
%% --------------------------------------------------------------------------
delta=0.2;

for l=6
   
    m=4;
    h=delta/m;
    n=1/h;
    
    %% A ------------------------------------------------------------------
    A=zeros(n,n);
    A(1,1)=2*h/3;
    A(1,n)=h/6;
    A(n,1)=h/6;
    
    for i=2:n
        A(i,i)=2*h/3;
    end
    
    for i=1:n-1
        A(i,i+1)=h/6;   
        A(i+1,i)=h/6;
    end
    

    
    
    %% B ------------------------------------------------------------------
    B=zeros(n+1,n);
    
   
    
    %band 1
    for i=1:n
        B(i,i)=h^2/6-m*h^2/2; %(x_i, x_i1)
    end
    
     %band 2
    for i=1:n
        B(i+1,i)=5*h^2/6-m*h^2/2;
    end
   
     %band 3--->m
     
     for j=3:m
         for i=j:n+1
             B(i,i-(j-1))=h^2;
         end
     end
   
        
    %band m+1
    for i=1:n-m+1
        B(m+i,i)=5*h^2/6;
    end   
     
    
    %band m+2
    
    for i=1:n-m
        B(m+i+1,i)=h^2/6;
    end
    
    % Add bands to 1:m rows and n+1 - th row
    
    %band 1
    for i=1:m+1
        B(i,n-m+i-1)=h^2/6;
    end
    
    
    
    %band 2
    for i=1:m
        B(i,n-m+i)=h^2*5/6;
    end
   
   
    %band 3-->m
    for i=1:m-2
         for j=1:m-i
             B(j,n-m+i+j)=h^2;
         end
     end
    
    

    
    %band m+1
    B(1,n)=5*h^2/6-m*h^2/2;
    
    %row n+1
    
    B(n+1,1)=h^2/6-m*h^2/2;
    
    
    % add coefficient
    
    B=-2/(delta^2)*B;   
    
  
    B(n+1,:)=[];
    
    
    %% RHS
   
    
   
     % Real solution
    ureal=@(x) sin(2*pi*x);  
    preal=@(x) -(2*((cos(2*pi*x) - cos(2*pi*(delta - x)))/(2*pi) + delta*sin(2*pi*x)))/delta^2;    
    f=@(x) (4*sin(2*pi*x)*(cos(2*pi*delta)/2 - delta^2*pi^2 + delta*pi*sin(2*pi*delta) - 1/2))/(delta^4*pi^2);
     
    F=Fgauss(h,f);
    c=mean(F);
    F=F-c;
        
    FF=[zeros(n,1);F];%2n
   
    
    MM=zeros(2*n,2*n);
 
    MM(1:n,1:n)=A;
    MM(n+1:2*n,1:n)=B';
    MM(1:n,n+1:2*n)=B;
    
    
    
    
    [xx,flag,relres] = minres(MM,FF); 
    %xx=MM\FF;
   
    %{
    c1=mean(xx(1:n));
    xx(1:n)=xx(1:n)-c1;
    
   
    c2=mean(xx(n+1:2*n));
    xx(1+n:2*n)=xx(1+n:2*n)-c2;
    
   %}
   

    figure
    
    plot(0:h:1,[xx(1:n);xx(1)],'y','LineWidth',2);
    hold on
    plot(0:h:1,preal(0:h:1),'k');
    legend('numerical-p', 'real-p')
    xlabel('x')
    ylabel('function value')
    title('NumSol-P')
    
    
    figure
    
    plot(h/2:h:1-h/2,xx(1+n:2*n),'y','LineWidth',2);
    hold on
    plot(h/2:h:1-h/2,ureal(h/2:h:1-h/2),'k');
    legend('numerical-u', 'real-u')
    xlabel('x')
    ylabel('function value')
    title('NumSol-U')
    
    
    %% err
    
    errP=getDeltaError(xx,f,delta,m);
   
    errU=getL2Error(xx,ureal,delta,m);
    
     
       
end


