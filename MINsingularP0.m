%function 
%input n delta
%output W

function [W]=MINsingularP0(n,delta)

h=1/n;


%%
%partI L2-norm
M0=eye(n,n);
M0=h*M0;

%partII delta-norm-reaction

H=[1, -1; -1, 1];

M1=zeros(n,n);

for i=1:n-1
    for j=i:i+1
        for k=i:i+1
            M1(j,k)=M1(j,k)+H(j-i+1,k-i+1);
        end
    end
end
M1(1,1)=M1(1,1)+1;
M1(1,n)=M1(1,n)-1;
M1(n,1)=M1(n,1)-1;
M1(n,n)=M1(n,n)+1;

M1=4/(3*delta)*M1;
%

M0=M0+M1;

%%
%N for L2-space
N=eye(n,n);
N=h*N;




%%

%SX & SY
SX=chol(M0);
%SX=SX';
SY=chol(N);
%SY=SY';






%%
%B
h1=-1*ones(n,1);
h12=ones(n-1,1);
B1=diag(h1);
B1=B1+diag(h12,-1);

B1(1,n)=B1(1,n)+1;


%Sigma=inv(SX)*B*inv(SY);
%Sigma=SX*B*SY;
%Sigma=SY*B'*SX;
%
Sigma=inv(SY)'*B1'*inv(SX);
s=svd(Sigma);
s_min=min(s);
W=s_min;

end


