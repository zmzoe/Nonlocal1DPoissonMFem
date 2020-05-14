function err=getError_sin_ph(Alpha,xx,f,delta,m)


h=delta/m;
n=1/h;

GaussP=[-0.7745967 0 0.7745967];                            %��˹��
GaussA=[0.5555556 0.8888889 0.5555556];                     %��˹ϵ��

x = 0:h:1;                                                  %����[0,1]
err=0;
F=zeros(length(x)-1,1);
for i=1:length(x)-1
    points = (x(i+1)-x(i))/2*GaussP + (x(i+1)+x(i))/2;      %����任 [-1,1]--->[x_i,x_i1] 
    
    [aa_i,bb_i,cc_i]=grad_psi(Alpha,delta,m,i);
   
    alpha_i=aa_i*xx(1:2*n-1);
    beta_i=bb_i*xx(1:2*n-1);
    gamma_i=cc_i*xx(1:2*n-1);
    
    alpha_x=@(x) alpha_i*x^2+beta_i*x+gamma_i;
    
    w=@(x) abs(f(x)-alpha_x(x));
    
    for k=1:3
        F(i) = F(i) + h/2*w(points(k))*GaussA(k);
    end
    
    err=err+F(i);
end
