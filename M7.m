%%  Main
% real u = sin (2*pi*x)
% real p = Du 
%      p =
%      2*pi*cos(2*pi*x)
% RHS  f = -D^2 u
%      f =
%      4*pi^2*sin(2*pi*x)
%% --------------------------------------------------------------------------

delta=0.1;
for l=2
    
    m=1;
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
    

    
    A
    %% B ------------------------------------------------------------------
    v1=-1*ones(n,1);
    B=diag(v1);
    
    v2=ones(n-1,1);
    B=B+diag(v2,-1);
     
    B(1,n)=1;
    
    
    
   
    
    
   % B=-B;

    
   
    
    %% RHS
    % Real solution
    ureal=@(x) sin(2*pi*x);
    
    preal=@(x) 2*pi*cos(2*pi*x);
    
    f=@(x) 4*pi^2*sin(2*pi*x);
     
    F=Fgauss(h,f);
    
   
   
    
    for i=1:n-1
        F(i)=F(i)-F(n);
    end
    F(n)=[];
    
        
    FF=[zeros(n-1,1);F];%2n-2
    
    
    E=eye(n);
    for i=1:n-1
        E(n,i)=-1;
        A=E'*A*E;
        B=E'*B*E;
        E=eye(n);      
    end
    
    MM=zeros(2*n-2,2*n-2);
 
    MM(1:n-1,1:n-1)=A(1:n-1,1:n-1);
    MM(n:2*n-2,1:n-1)=-B(1:n-1,1:n-1)';
    MM(1:n-1,n:2*n-2)=B(1:n-1,1:n-1);
    
    xx=MM\FF;

    
    figure 
    plot(xx);
  

    figure
    
    plot(0:h:1,[xx(1:n-1);-sum(xx(1:n-1));xx(1)],'y','LineWidth',2);
    hold on
    plot(0:h:1,preal(0:h:1),'k');
   
    
    
    
    figure
    
    plot(h/2:h:1-h/2,[xx(n:2*n-2);-sum(xx(n:2*n-2))],'y','LineWidth',2);
    hold on
    plot(h/2:h:1-h/2,ureal(h/2:h:1-h/2),'k');
   
    
    
    %% err
    
   
     
       
end


