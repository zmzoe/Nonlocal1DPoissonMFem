function G=GgaussDirectPolyHH(h,g,delta)


GaussP=[-0.9061798 -0.5384693 0 0.5384693 0.9061798];                            %高斯点 9次代数精度
GaussA=[0.2369269 0.4786287 0.5688889 0.4786287 0.236926];                       %高斯系数

XX = 0:h:1;                                                  %区间[0,1]
n = length(XX)-1;
G=zeros(2*n,1);

%% odd  被积函数 g*(phi_2i-1 - phi_2i)
for i=1:n
    points = h/2*GaussP + (XX(i+1)+XX(i))/2;                  %区间变换(寻找[x_i,x_i1]上相应的积分点)
    
   for k=1:5
       if points(k)<1/2-delta && points(k)>=0
           gg=@(x) g{1}(x);
       elseif points(k)<1/2 && points(k)>=1/2-delta
               gg=@(x) g{2}(x);
       elseif points(k)<1-delta && points(k)>=1/2
               gg=@(x) g{3}(x);
       elseif points(k)<=1 && points(k)>=1-delta
               gg=@(x) g{4}(x);
                   
       end
                    
    w =@(x) gg(x)*(XX(i+1)-x)/h;                             %被积函数 g*(phi_2i-1)
    
    G(2*i-1) = G(2*i-1) + h/2*w(points(k))*GaussA(k);
    
   end
end



%% even  被积函数 g*( phi_2i)
for i=1:n
    points = h/2*GaussP + (XX(i+1)+XX(i))/2;                  %区间变换(寻找[x_i,x_i1]上相应的积分点)
    
    
                            
     
    for k=1:5
       if points(k)<1/2-delta && points(k)>=0
           gg=@(x) g{1}(x);
       elseif points(k)<1/2 && points(k)>=1/2-delta
               gg=@(x) g{2}(x);
       elseif points(k)<1-delta && points(k)>=1/2
                   gg=@(x) g{3}(x);
       elseif points(k)<=1 && points(k)>=1-delta
                   gg=@(x) g{4}(x);
                   
       end
                    
        w =@(x) gg(x)*(x-XX(i))/h;         %被积函数 g*(phi_2i)
        G(2*i) = G(2*i) + h/2*w(points(k))*GaussA(k);
    end
end



