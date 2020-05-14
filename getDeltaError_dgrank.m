
%   For the basis with both compatibility and boundary condition
%   Input : 
%   1) The basis (Alpha is the transformation matrix
%               between any psi and the half hat function basis phi)
%   2) xx (numerical results for both p and u
%        n-1 + n-1 or 2*n-1 + n-1
%
%   3) f (notice that Gp=f)
%
%   4) delta,m 

%   Output :
%   error \|G_delta^*(p-p_h)\|_0

%------------------------------------------------------
%   How to realize this? 
%   use gauss integration formula, first of all find gauss points in [-1,1]
%   and corresponding coefficients. Second, our domain \Omega=[0,1] can be
%   seprated to n intervals I_i i=1:n  so we want to find the gauss point
%   at every I_i for every i.
%   
%   Next,grad_psi is to compute G_delta Psi over I_i for each i
%   
%
%--------------------------------------------------------
function err=getDeltaError_dgrank(Alpha,xx,f,delta,m)


h=delta/m;
n=1/h;

GaussP=[-0.7745967 0 0.7745967];                            %高斯点
GaussA=[0.5555556 0.8888889 0.5555556];                     %高斯系数

x = 0:h:1;                                                  %区间[0,1]
err=0;
F=zeros(length(x)-1,1);
for i=1:length(x)-1
    points = (x(i+1)-x(i))/2*GaussP + (x(i+1)+x(i))/2;      %区间变换 [-1,1]--->[x_i,x_i1] 
    
    [aa_i,bb_i,cc_i]=grad_psi(Alpha,delta,m,i);
   
    alpha_i=aa_i*xx(1:2*n-1);
    beta_i=bb_i*xx(1:2*n-1);
    gamma_i=cc_i*xx(1:2*n-1);
    
    alpha_x=@(x) alpha_i*x^2+beta_i*x+gamma_i;
    
    w=@(x) (f(x)+alpha_x(x))^2;
    
    for k=1:3
        F(i) = F(i) + h/2*w(points(k))*GaussA(k);
    end
    
    err=err+F(i);
end
