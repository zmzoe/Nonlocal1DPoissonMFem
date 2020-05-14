%function [W]=M1d5singular(delta)
delta=0.5;
h=delta;
d=delta;
%M=1 h=delta/M delta=0.5
%% ~~~~~~~~~computation of \|p\|_delta^2~~~~~~~~~~~~~~~~~~~~~
H11=(2*h/3+40*h/(3*d^2)-8/d+2/h-8*h^2/(d^3)+(8*h^3)/(5*d^4))*ones(2,1);
M=diag(H11);

M(2,1)=h/3+8/d-2/h+(8*h^2)/(d^3)-(8*h^3)/(5*d^4)-(40*h)/(3*d^2);
M(1,2)=M(2,1);




%% ~~~~~~~~~ the computation of \|p\|_0^2~~~~~~~~~~~~~~~~~~~~~
          
%N for L2-space
N=eye(2,2);
N=h*N;




%%
%SX & SY
SX=chol(M);
%SX=SX';
SY=chol(N);
%SY=SY';


%%
%B
B=(1+(2*h^2)/(3*d^2)-(2*h)/d)*ones(2,2);
B(2,1)=-B(2,1);
B(1,2)=-B(1,2);


Sigma=inv(SY)'*B'*inv(SX);
s=svd(Sigma);
s_min=min(s);
W=s_min;
W


