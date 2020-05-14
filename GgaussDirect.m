

function G=GgaussDirect(h,g)

XX = 0:h:1;                                                  %区间[0,1]
n = length(XX)-1;
G=zeros(2*n,1);

%% odd  被积函数 g*(phi_2i-1 - phi_2i)
for i=1:n
   

    wodd =@(y) g(y).*(XX(i+1)-y)/h;                              %被积函数 g*(phi_2i-1)
    
    weven =@(y) g(y).*(y-XX(i))/h;   
     
    G(2*i-1) = integral(@(y)wodd(y), XX(i), XX(i+1),'AbsTol', 1e-12); 
    G(2*i) = integral(@(y)weven(y), XX(i), XX(i+1),'AbsTol', 1e-12); 
    
    
end







