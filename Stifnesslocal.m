% Frame
% Import alpha_i, beta_i .   row vector
% Output  (Dalpha,beta)

function [DAlphaBeta]=Stifnesslocal(alpha,beta, n )
h=1/n;

A=zeros(2*n,2*n);

for i=1:n
    A(2*i-1,2*i-1)=1/h;
    A(2*i-1,2*i)=-1/h;
    
    A(2*i,2*i-1)=-1/h;
    A(2*i,2*i)=1/h;
end



DAlphaBeta = alpha*A*beta';
end