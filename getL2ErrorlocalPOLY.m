function err=getL2ErrorlocalPOLY(uu,u,h)


GaussP=[-0.9061798 -0.5384693 0 0.5384693 0.9061798];                            %高斯点 9次代数精度
GaussA=[0.2369269 0.4786287 0.5688889 0.4786287 0.236926];                       %高斯系数


%GaussP=[-0.7745967 0 0.7745967];                            %高斯点
%GaussA=[0.5555556 0.8888889 0.5555556];                     %高斯系数

x = 0:h:1;                                                  %区间[0,1]
err=0;
F=zeros(length(x)-1,1);
for i=1:length(x)-1
    points = h/2*GaussP + (x(i+1)+x(i))/2;                  %区间变换
    
    for k=1:5
       if points(k)<1/2 && points(k)>=0
               ureal=@(x) u{1}(x);
       elseif points(k)<=1 && points(k)>=1/2
               ureal=@(x) u{2}(x);
               
       end
    
    w=@(x) (ureal(x)-uu(i))^2;
    
   
        F(i) = F(i) + h/2*w(points(k))*GaussA(k);
    end
    
    err=err+F(i);
end