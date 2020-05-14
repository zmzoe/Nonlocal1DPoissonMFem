% Frame
% Import alpha_i, beta_i .   row vector
% Output  (Dalpha,beta)

function [DAlphaBeta]=Ponelocal(alpha,beta, n )


B=zeros(2*n,2*n);

for i=1:n
    B(2*i-1,2*i-1)=-1/2;
    B(2*i-1,2*i)=-1/2;
    
    B(2*i,2*i-1)=1/2;
    B(2*i,2*i)=1/2;
end



DAlphaBeta = alpha*B*beta';
end