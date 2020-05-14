
   l=2
    m=2^l;
    h=delta/m;
    n=1/h;
    
   
    
    

    
    
    %% B ------------------------------------------------------------------
    B=zeros(n,n);
    
   
    
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
    
    
    
    
    % add coefficient
    
    B=-2/(delta^2)*B;   
    
  
  
    