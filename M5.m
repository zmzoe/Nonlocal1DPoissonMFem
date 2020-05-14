
%%  Main
% real u = sin (pi*x)
% real p = Gu 
%      p = 2/(delta^2)*1/pi*(cos(pi*(x-delta))-cos(pi*x))-2/(delta)*sin(pi*x)
% RHS  f = -G^*G u
%      G^*G u = 4*sin(pi*x)/(delta^2)+16*(sin((pi/2)*delta))^2*sin(pi*x)/(delta^4*pi^2)-8*sin(pi*delta)*sin(pi*x)/(delta^3*pi)
%
%
%% --------------------------------------------------------------------------
delta=0.2;

for l=6
    m=2^l;
    h=delta/m;
    n=1/h;
    
    %% A ------------------------------------------------------------------
    A=zeros(n+1,n+1);
    A(1,1)=h/3;
    A(n+1,n+1)=h/3;
    
    for i=2:n
        A(i,i)=2*h/3;
    end
    
    for i=1:n
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
    
  
    
    
    
    %% RHS
    % Real solution
    ureal=@(x) sin(pi*x);
    
    preal=@(x) 2/(delta^2)*1/pi*(cos(pi*(x-delta))-cos(pi*x))-2/(delta)*sin(pi*x);
    
    f=@(x) -(4*sin(pi*x)/(delta^2)+16*(sin((pi/2)*delta))^2*sin(pi*x)/(delta^4*pi^2)-8*sin(pi*delta)*sin(pi*x)/(delta^3*pi));   
     
    F=Fgauss(h,f);
    
    % f without CC
    
    FF=[zeros(n+1,1);F;0;0;0];
    
    
  
    
    %% form the whole matrix
    
    MM=zeros(2*n+1,2*n+1);
 
    MM(1:n+1,1:n+1)=A;
    MM(n+2:2*n+1,1:n+1)=B';
    MM(1:n+1,n+2:2*n+1)=B;
    
    
    
    ZZ1=[ones(1,n),zeros(1,n+1)];
    ZZ2=[1,zeros(1,n-1),-1,zeros(1,n)];
    ZZ3=[zeros(1,n+1),ones(1,n)];
    NN=[MM,ZZ1',ZZ2',ZZ3';ZZ1,0,0,0;ZZ2,0,0,0;ZZ3,0,0,0];
   
    
    %% solve
    xx=NN\FF;
    node=0:h:1;
    
    figure
    plot(h/2:h:1-h/2,xx(n+2:2*n+1));
    hold on 
    plot(h/2:h:1-h/2,ureal(h/2:h:1-h/2));
    
    figure
    plot(0:h:1,xx(1:n+1));
    hold on 
    plot(0:h:1,preal(0:h:1));%}
    
    %% ERROR
    
  

    
    
    
    
   
   
     
       
end


