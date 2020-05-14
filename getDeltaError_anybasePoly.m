function err=getDeltaError_anybasePoly(Alpha,xx,f,delta,m)


h=delta/m;
n=1/h;

%GaussP=[-0.7745967 0 0.7745967];                            %��˹�� 5�δ�������
%GaussA=[0.5555556 0.8888889 0.5555556];                     %��˹ϵ��


GaussP=[-0.9061798 -0.5384693 0 0.5384693 0.9061798];                            %��˹�� 9�δ�������
GaussA=[0.2369269 0.4786287 0.5688889 0.4786287 0.236926];                       %��˹ϵ��
x = 0:h:1;                                                  %����[0,1]
err=0;
F=zeros(length(x)-1,1);
for i=1:length(x)-1
    points = (x(i+1)-x(i))/2*GaussP + (x(i+1)+x(i))/2;      %����任 [-1,1]--->[x_i,x_i1] 
    
    [aa_i,bb_i,cc_i]=grad_psi(Alpha,delta,m,i);
   
    alpha_i=aa_i*xx(1:n-1);
    beta_i=bb_i*xx(1:n-1);
    gamma_i=cc_i*xx(1:n-1);
    
    alpha_x=@(x) alpha_i*x^2+beta_i*x+gamma_i;
    
    for k=1:5
       if points(k)<delta && points(k)>=0
               ff=@(x) f{1}(x);
       elseif points(k)< 1/2-delta && points(k)>=delta
               ff=@(x) f{2}(x);
       elseif points(k)< 1/2 && points(k)>=1/2-delta
               ff=@(x) f{3}(x);
       elseif points(k)<1/2+delta && points(k)>=1/2
               ff=@(x) f{4}(x);
       elseif points(k)<1-delta && points(k)>=1/2+delta
               ff=@(x) f{5}(x);
       elseif points(k)<=1 && points(k)>=1-delta
               ff=@(x) f{6}(x);
                   
       end
       
       w=@(x) (ff(x)+alpha_x(x))^2;
       
        F(i) = F(i) + h/2*w(points(k))*GaussA(k);
    end
    
    err=err+F(i);
end
