function G=Ggauss(h,g)
GaussP=[-0.7745967 0 0.7745967];                            %高斯点
GaussA=[0.5555556 0.8888889 0.5555556];                     %高斯系数

XX = 0:h:1;                                                  %区间[0,1]
n = length(XX)-1;
G=zeros(2*n-1,1);

%% odd  被积函数 g*(phi_2i-1 - phi_2i)
for i=1:n
    points = h/2*GaussP + (XX(i+1)+XX(i))/2;                  %区间变换(寻找[x_i,x_i1]上相应的积分点)
    
    
    w =@(x) g(x)*(XX(i+1)+XX(i)-2*x)/h;                              %被积函数 g*(phi_2i-1 - phi_2i)
     
    for k=1:3
        G(2*i-1) = G(2*i-1) + h/2*w(points(k))*GaussA(k);
    end
end

%% even 被积函数 g*(phi_2i)在[x_i,x_i1], 减去 g*(phi_2i1)在[x_i1,x_i2]

% part 1-------被积函数 g*(phi_2i)在[x_i,x_i1]

G1 = zeros(n-1,1);
G2 = zeros(n-1,1);
for i=1:n-1
    points = h/2*GaussP + (XX(i+1)+XX(i))/2;                  %区间变换(寻找[x_i,x_i1]上相应的积分点)
    
   
    w =@(x) g(x)*(x-XX(i))/h;                                       % 被积函数 g*(phi_2i)
     
    for k=1:3
        G1(i) = G1(i) + h/2*w(points(k))*GaussA(k);
    end
end

% part 2-------被积函数 g*(phi_2i1)在[x_i1,x_i2]
for i=1:n-1
    points = h/2*GaussP + (XX(i+2)+XX(i+1))/2;                  %区间变换(寻找[x_i,x_i1]上相应的积分点)
    
   
    w =@(x) g(x)*(XX(i+2)-x)/h;                                       %被积函数 g*(phi_2i1)
    
    for k=1:3
        G2(i) = G2(i) + h/2*w(points(k))*GaussA(k);
    end
end

for i=1:n-1
    G(2*i)= G1(i)-G2(i);
end

