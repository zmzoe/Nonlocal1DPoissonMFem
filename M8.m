 %Main
% real u = sin (2*pi*x)
% real p = Gu 
%      p =-(2*((cos(2*pi*x) - cos(2*pi*(delta - x)))/(2*pi) + delta*sin(2*pi*x)))/delta^2
%      
% RHS  f = -G^*G u
%      f =(4*sin(2*pi*x)*(cos(2*pi*delta)/2 - delta^2*pi^2 + delta*pi*sin(2*pi*delta) - 1/2))/(delta^4*pi^2)
%      
%
%% --------------------------------------------------------------------------
delta=0.2;

for l=4
   
    m=2^l;
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
    
    for i=1:n-1
        F(i)=F(i)-F(n);
    end
    F(n)=[];
    
        
    FF=[zeros(n-1,1);F];%2n-2
   
    MM=zeros(2*n,2*n);
 
    MM(1:n,1:n)=A;
    MM(n+1:2*n,1:n)=B';
    MM(1:n,n+1:2*n)=B;
    
   
    %% boundary condition and compatibility condition 
    
    E1=eye(2*n);
    E2=eye(2*n);
    for i=1:n-1
        E1(n,i)=-1;
        MM=E1'*MM*E1;
        E2(2*n,i+n)=-1;
        MM=E2'*MM*E2;
        E1=eye(2*n); 
        E2=eye(2*n); 
    end
    
    
    MM(n,:)=[];
    MM(:,n)=[];
    
    MM(2*n-1,:)=[];
    MM(:,2*n-1)=[];
    
       %xx=MM\FF;
      [xx,flag,relres] = gmres(MM,FF);   
    
    
      xx=[xx(1:n-1);-sum(xx(1:n-1));xx(1);xx(n:2*n-2);-sum(xx(n:2*n-2))];
    
    
    
    
    
    

   
    
     figure
    
    plot(0:h:1,[xx(1:n-1);-sum(xx(1:n-1));xx(1)],'y','LineWidth',2);
    hold on
    plot(0:h:1,preal(0:h:1),'k');
   
    
    
    
    figure
    
    plot(h/2:h:1-h/2,[xx(n:2*n-2);-sum(xx(n:2*n-2))],'y','LineWidth',2);
    hold on
    plot(h/2:h:1-h/2,ureal(h/2:h:1-h/2),'k');
   

   
    %% err
    
    errP=getDeltaError(xx,f,delta,m);
   
    errU=getL2Error(xx,ureal,delta,m);
    
    
   
   
    
    
    
     
       
end

