

%(x_i-h/2,x_i+h/2)

function F=FgaussMove(h,f)
GaussP=[-0.7745967 0 0.7745967];                            %高斯点
GaussA=[0.5555556 0.8888889 0.5555556];                     %高斯系数

x = 0:h:1;                                                  %区间[0,1]
F=zeros(length(x)-1,1);
for i=1:length(x)-1
    points = h/2*GaussP + (x(i)+x(i))/2;                  %区间变换 (b-a)/2*point+(a+b)/2
    
    for k=1:3
        F(i) = F(i) + h/2*f(points(k))*GaussA(k);
    end
end
