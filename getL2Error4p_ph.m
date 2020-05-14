function err=getL2Error4p_ph(p_h,p,h)

XX= 0:h:1; 
n=length(XX)-1;
err=0;
F=zeros(n,1);
for i=1:n
    p_K=@(y) p_h(i).*(XX(i+1)-y)/h+p_h(i+1).*(y-XX(i))/h;
    w=@(y) (p(y)-p_K(y)).^2;
    
    F(i) = integral(@(y)w(y), XX(i), XX(i+1),'AbsTol', 1e-12);
    
    err=err+F(i);
end