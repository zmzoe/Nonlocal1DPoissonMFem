function err=getDeltaError_anybase(Alpha,xx,HH,delta,m)


h=delta/m;
n=1/h;

%GaussP=[-0.7745967 0 0.7745967];                            %高斯点 5次代数精度
%GaussA=[0.5555556 0.8888889 0.5555556];                     %高斯系数


GaussP=[-0.9061798 -0.5384693 0 0.5384693 0.9061798];                            %高斯点 9次代数精度
GaussA=[0.2369269 0.4786287 0.5688889 0.4786287 0.236926];                       %高斯系数
x = 0:h:1;                                                  %区间[0,1]
err=0;
F=zeros(length(x)-1,1);
for i=1:length(x)-1
    points = (x(i+1)-x(i))/2*GaussP + (x(i+1)+x(i))/2;      %区间变换 [-1,1]--->[x_i,x_i1] 
    
    [aa_i,bb_i,cc_i]=grad_psi(Alpha,delta,m,i);
   
    alpha_i=aa_i*xx(1:n-1);
    beta_i=bb_i*xx(1:n-1);
    gamma_i=cc_i*xx(1:n-1);
    
    alpha_x=@(x) alpha_i*x^2+beta_i*x+gamma_i;
    
    w=@(x) (-HH(x)-alpha_x(x))^2;
    
    for k=1:5
        F(i) = F(i) + h/2*w(points(k))*GaussA(k);
    end
    
    err=err+F(i);
end
