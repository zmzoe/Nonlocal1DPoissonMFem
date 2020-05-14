function err=getL2Error_local_p(pp,p,h)


GaussP=[-0.9061798 -0.5384693 0 0.5384693 0.9061798];                            %��˹�� 9�δ�������
GaussA=[0.2369269 0.4786287 0.5688889 0.4786287 0.236926];                       %��˹ϵ��

XX= 0:h:1;                                                  %����[0,1]
err=0;
F=zeros(length(XX)-1,1);% n
for i=1:length(XX)-1
    points = h/2*GaussP + (XX(i+1)+XX(i))/2;                  %����任
    
    
    p_K=@(y) pp(i)*(XX(i+1)-y)/h+pp(i+1)*(y-XX(i))/h;
    w=@(y) (p(y)-p_K(y))^2;
    
    for k=1:5
        F(i) = F(i) + h/2*w(points(k))*GaussA(k);
    end
    
    err=err+F(i);
end