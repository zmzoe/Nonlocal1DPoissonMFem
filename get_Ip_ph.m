function err=get_Ip_ph(Alpha,p_h,p_r,delta,m)


h  = delta/m;
n  = 1/h;     
XX = 0:h:1;
err= 0;
F  = zeros(n,1);

for i=1:n
    [aa_i,bb_i,cc_i]=grad_psi(Alpha,delta,m,i);
   
    alpha1 = aa_i*p_h(1:n);
    beta1= bb_i*p_h(1:n);
    gamma1=cc_i*p_h(1:n);

    alpha_p_h=@(y) alpha1*y.^2+beta1*y+gamma1;
    
    alpha2 = aa_i*p_r(1:n);
    beta2= bb_i*p_r(1:n);
    gamma2=cc_i*p_r(1:n);
    alpha_p_r=@(y) alpha2*y.^2+beta2*y+gamma2;
    
    w=@(y) (alpha_p_r(y)-alpha_p_h(y)).^2;

    F(i) = integral(@(y)w(y), XX(i), XX(i+1),'AbsTol', 1e-12);
    
    err=err+F(i);
end
