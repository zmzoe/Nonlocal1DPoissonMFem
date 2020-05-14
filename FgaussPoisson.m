function F=FgaussPoisson(h,f)


GaussP=[-0.9061798 -0.5384693 0 0.5384693 0.9061798];                            %高斯点 9次代数精度
GaussA=[0.2369269 0.4786287 0.5688889 0.4786287 0.236926];   



XX = 0:h:1;                                                  %区间[0,1]

F=zeros(length(XX)-1,1);
n=length(XX)-1;
%% Part I
for i=1:length(XX)-1
    points = h/2*GaussP + (XX(i+1)+XX(i))/2;                  %区间变换
    
    phi_minus=@(y) (XX(i+1)-y)/h;
    fphi_minus=@(y) f(y)*phi_minus(y);
    
    for k=1:5
        F(i) = F(i) + h/2*fphi_minus(points(k))*GaussA(k);
    end
end
%% Part II

for i=2:length(XX)-1
    points = h/2*GaussP + (XX(i)+XX(i-1))/2;                  %区间变换
    
    phi_plus=@(y) (y-XX(i))/h;
    fphi_plus=@(y) f(y)*phi_plus(y);
    for k=1:5
        F(i) = F(i) + h/2*fphi_plus(points(k))*GaussA(k);
    end
end

%%
points = h/2*GaussP + (XX(n)+XX(n+1))/2;
phi1plus=@(y) (y-XX(n))/h;
fphi1plus=@(y) f(y)*phi1plus(y);
for k=1:5
        F(1) = F(1) + h/2*fphi1plus(points(k))*GaussA(k);
end


