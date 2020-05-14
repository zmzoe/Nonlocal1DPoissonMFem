%
% Frame
% Import alpha_i, beta_i .   row vector
% Output (alphi,beta)

function [InnerP]=innerproductlocal(alpha,beta,n)
%%  Parameter
h=1/n;

vv = h/3*ones(2*n,1);

uu = zeros(2*n-1,1);

for i = 1:n

    uu(2*i-1,1) = h/6; % forget 0

end


D  = diag(vv);
D  = D+diag(uu,-1)+diag(uu,1);


    
    
    

InnerP  = alpha*D*beta';
end