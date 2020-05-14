function err=getL2Errorlocal(u_h,u,h)


XX = 0:h:1;                                                  
n=length(XX)-1;
err=0;
F=zeros(n,1);
for i=1:n
    
    w=@(y) (u(y)-u_h(i)).^2;
    
    F(i) = integral(@(y)w(y), XX(i), XX(i+1),'AbsTol', 1e-12);
    
    err=err+F(i);
end