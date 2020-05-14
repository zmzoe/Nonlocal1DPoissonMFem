function err=getDeltaError_anybase1(Alpha,xx,HH,delta,m)


h=delta/m;
n=1/h;

%GaussP=[-0.7745967 0 0.7745967];                            %��˹�� 5�δ�������
%GaussA=[0.5555556 0.8888889 0.5555556];                     %��˹ϵ��


GaussP=[-0.9061798 -0.5384693 0 0.5384693 0.9061798];                            %��˹�� 9�δ�������
GaussA=[0.2369269 0.4786287 0.5688889 0.4786287 0.236926];                       %��˹ϵ��
XX = 0:h:1;                                                  %����[0,1]
err=0;
F=zeros(n,1);
for i=1:n
    points = (XX(i+1)-XX(i))/2*GaussP + (XX(i+1)+XX(i))/2;      %����任 [-1,1]--->[x_i,x_i1] 
    
    [aa_i,bb_i,cc_i]=grad_psi(Alpha,delta,m,i);
   
    alpha_i=aa_i*xx(1:n);
    beta_i=bb_i*xx(1:n);
    gamma_i=cc_i*xx(1:n);
    
    alpha_x=@(y) alpha_i*y.^2+beta_i*y+gamma_i;
    
    w=@(y) (-HH(y)-alpha_x(y)).^2;
    
    for k=1:5
        F(i) = F(i) + h/2*w(points(k))*GaussA(k);
    end
    
    err=err+F(i);
end
