function err=getL2Error(uu,u,delta,m)


h=delta/m;
n=1/h;

GaussP=[-0.7745967 0 0.7745967];                            %高斯点
GaussA=[0.5555556 0.8888889 0.5555556];                     %高斯系数

x = 0:h:1;                                                  %区间[0,1]
err=0;
F=zeros(length(x)-1,1);
for i=1:length(x)-1
    points = h/2*GaussP + (x(i+1)+x(i))/2;                  %区间变换
    
    
    
    w=@(x) (u(x)-uu(i))^2;
    
    for k=1:3
        F(i) = F(i) + h/2*w(points(k))*GaussA(k);
    end
    
    err=err+F(i);
end