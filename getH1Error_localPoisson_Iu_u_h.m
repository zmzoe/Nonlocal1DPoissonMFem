function err=getH1Error_localPoisson_Iu_u_h(pp,ureal,h)



GaussP=[-0.9061798 -0.5384693 0 0.5384693 0.9061798];                            %高斯点 9次代数精度
GaussA=[0.2369269 0.4786287 0.5688889 0.4786287 0.236926];   

x= 0:h:1;                                                  %区间[0,1]
err=0;
F=zeros(length(x)-1,1);% n
for i=1:length(x)-1
    points = h/2*GaussP + (x(i+1)+x(i))/2;                  %区间变换
    
    
    p_K=@(y) (pp(i+1)-pp(i))/h;
    Iu=@(y) (ureal(i*h)-ureal((i-1)*h))/h;
    w=@(y) (Iu(y)-p_K(y))^2;
    
    for k=1:5
        F(i) = F(i) + h/2*w(points(k))*GaussA(k);
    end
    
    err=err+F(i);
end