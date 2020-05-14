function err=getL2Error_local_pPOLY(pp,p,h,delta)


%GaussP=[-0.7745967 0 0.7745967];                            %��˹��
%GaussA=[0.5555556 0.8888889 0.5555556];                     %��˹ϵ��

GaussP=[-0.9061798 -0.5384693 0 0.5384693 0.9061798];                            %��˹�� 9�δ�������
GaussA=[0.2369269 0.4786287 0.5688889 0.4786287 0.236926];                       %��˹ϵ��

x= 0:h:1;                                                  %����[0,1]
err=0;
F=zeros(length(x)-1,1);% n
for i=1:length(x)-1
    points = h/2*GaussP + (x(i+1)+x(i))/2;                  %����任
    
    for k=1:5
       if points(k)<delta && points(k)>=0
               preal=@(y) p{1}(y);
       elseif points(k)<1/2 && points(k)>=delta
               preal=@(y) p{2}(y);
       elseif points(k)<1/2+delta && points(k)>=1/2
               preal=@(y) p{3}(y);
       elseif points(k)<=1 && points(k)>=1/2+delta
               preal=@(y) p{4}(y);
                   
       end
    p_K=@(y) pp(i)*(x(i+1)-y)/h+pp(i+1)*(y-x(i))/h;
    w=@(y) (preal(y)-p_K(y))^2;
    
        F(i) = F(i) + h/2*w(points(k))*GaussA(k);
    end
    
    err=err+F(i);
end