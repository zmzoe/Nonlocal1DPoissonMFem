function err=getL2ErrorInter(uh,ureal,h)


GaussP=[-0.9061798 -0.5384693 0 0.5384693 0.9061798];                            %高斯点 9次代数精度
GaussA=[0.2369269 0.4786287 0.5688889 0.4786287 0.236926];                       %高斯系数

XX= 0:h:1;                                                       %区间[0,1]
err=0;
F=zeros(length(XX)-1,1);% n
for i=1:length(XX)-1
    points = h/2*GaussP + (XX(i+1)+XX(i))/2;                     %区间变换
    
    
    p_K=@(y) uh(i)*(XX(i+1)-y)/h+uh(i+1)*(y-XX(i))/h;             %数值解
    p_h=@(y) ureal((i-1)*h)*(XX(i+1)-y)/h+ureal(i*h)*(y-XX(i))/h; %插值
    w=@(y) (p_h(y)-p_K(y))^2;
    
    for k=1:5
        F(i) = F(i) + h/2*w(points(k))*GaussA(k);
    end
    
    err=err+F(i);
end