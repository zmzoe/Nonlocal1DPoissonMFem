% jump 

% for any function yita, compute [yita]_x_i1

%alpha=1*2n
function [Jumpoverxi1,vecD]=jump(alpha, delta, m, i)
%%  Parameter

h = delta/m;
n = 1/h;

%% 

vecD=zeros(2*n,1);

    vecD(2*i,1)=-1;
    vecD(2*i+2,1)=1;


Jumpoverxi1=alpha*vecD;
