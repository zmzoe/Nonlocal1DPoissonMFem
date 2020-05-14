function err=getDeltaError4anybase(Alpha,p_h,f,delta,m)


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

    alpha_p=@(y) alpha1*y.^2+beta1*y+gamma1;
    
    w=@(y) (-f(y)-alpha_p(y)).^2;

    F(i) = integral(@(y)w(y), XX(i), XX(i+1),'AbsTol', 1e-12);
    
    err=err+F(i);
end
