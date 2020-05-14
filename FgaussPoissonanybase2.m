function F=FgaussPoissonanybase2(h,f)

XX = 0:h:1;                                                  %Çø¼ä[0,1]
n=length(XX)-1;

F=zeros(2*n,1);

for i=1:n
   
    phi_even=@(y) (XX(i+1)-y)/h;
    fphi_even=@(y) f(y).*phi_even(y);
    
    phi_odd=@(y) (y-XX(i))/h;
    fphi_odd=@(y) f(y).*phi_odd(y);
    
   
        F(2*i-1) = integral(@(y)fphi_even(y), (i-1)*h, i*h, ...
        'AbsTol', 1e-12);
        F(2*i) = integral(@(y)fphi_odd(y), (i-1)*h, i*h, ...
        'AbsTol', 1e-12);

end
