%
% Frame
% Import alpha_i, beta_i .   row vector
% Output (d alphi/dx,d beta/dx)

function [InnerGP]=innerDPproductlocal(alpha,beta,n)
%%  Parameter
h=1/n;

D=zeros(2*n,2*n);

for i = 1:n

    D(2*i-1,2*i-1) =  1/h;
    D(2*i-1,2*i)   = -1/h;
    D(2*i,2*i-1)   = -1/h;
    D(2*i,2*i)     =  1/h;
    

end


    
    
    

InnerGP  = alpha*D*beta';
end