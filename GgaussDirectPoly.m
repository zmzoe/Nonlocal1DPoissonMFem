function G=GgaussDirectPoly(h,g,delta)


GaussP=[-0.9061798 -0.5384693 0 0.5384693 0.9061798];                            %��˹�� 9�δ�������
GaussA=[0.2369269 0.4786287 0.5688889 0.4786287 0.236926];                       %��˹ϵ��

XX = 0:h:1;                                                  %����[0,1]
n = length(XX)-1;
G=zeros(2*n,1);

%% odd  �������� g*(phi_2i-1 - phi_2i)
for i=1:n
    points = h/2*GaussP + (XX(i+1)+XX(i))/2;                  %����任(Ѱ��[x_i,x_i1]����Ӧ�Ļ��ֵ�)
    
   for k=1:5
       if points(k)<delta && points(k)>=0
               gg=@(x) g{1}(x);
       elseif points(k)< 1/2-delta && points(k)>=delta
               gg=@(x) g{2}(x);
       elseif points(k)< 1/2 && points(k)>=1/2-delta
               gg=@(x) g{3}(x);
       elseif points(k)<1/2+delta && points(k)>=1/2
               gg=@(x) g{4}(x);
       elseif points(k)<1-delta && points(k)>=1/2+delta
               gg=@(x) g{5}(x);
       elseif points(k)<=1 && points(k)>=1-delta
               gg=@(x) g{6}(x);
                   
       end
                    
    w =@(x) gg(x)*(XX(i+1)-x)/h;                             %�������� g*(phi_2i-1)
    
    G(2*i-1) = G(2*i-1) + h/2*w(points(k))*GaussA(k);
    
   end
end



%% even  �������� g*( phi_2i)
for i=1:n
    points = h/2*GaussP + (XX(i+1)+XX(i))/2;                  %����任(Ѱ��[x_i,x_i1]����Ӧ�Ļ��ֵ�)
    
    
                            
     
    for k=1:5
       if points(k)<delta && points(k)>=0
               gg=@(x) g{1}(x);
       elseif points(k)< 1/2-delta && points(k)>=delta
               gg=@(x) g{2}(x);
       elseif points(k)< 1/2 && points(k)>=1/2-delta
               gg=@(x) g{3}(x);
       elseif points(k)<1/2+delta && points(k)>=1/2
               gg=@(x) g{4}(x);
       elseif points(k)<1-delta && points(k)>=1/2+delta
               gg=@(x) g{5}(x);
       elseif points(k)<=1 && points(k)>=1-delta
               gg=@(x) g{6}(x);
                   
       end
                    
        w =@(x) gg(x)*(x-XX(i))/h;         %�������� g*(phi_2i)
        G(2*i) = G(2*i) + h/2*w(points(k))*GaussA(k);
    end
end



