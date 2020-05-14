%function 
%input n delta
%output W

%function [W]=P1singular(n,delta)
delta=0.001;
n=5;
h=1/n;
%delta=1/4*h;%change
%basic information & variables

%%
%M for delta-space
%Part 1 delta-norm-nonlocal
H=zeros(3,3);
H(1,1:3)=[(15*h-7*delta)/(15*h^2), (3*delta-5*h)/(5*h^2), -(2*delta)/(15*h^2)];
H(2,2:3)=[(15*h-8*delta)/(15*h^2), -delta/(15*h^2)];
H(3,3)=delta/(5*h^2);

C=triu(H,1);
C=C';

H=H+C;

M=zeros(n+1,n+1);

for i=1:n-1
    for j=i:i+2
        for k=i:i+2
            M(j,k)=M(j,k)+H(j-i+1,k-i+1);
        end
    end
end
%% ~~~~~~~~~~~~~please consider boundary conditions~~~~~~~~~~~~~~~~~~~~~

M(1,1)=M(1,1)+H(3,3);

M(1,n)=M(1,n)+H(1,3);
M(1,n+1)=M(1,n+1)+H(2,3);

M(n,1)=M(n,1)+H(3,1);
M(n+1,1)=M(n+1,1)+H(3,2);

M(n,n)=M(n,n)+H(1,1);
M(n,n+1)=M(n,n+1)+H(1,2);
M(n+1,n)=M(n+1,n)+H(2,1);
M(n+1,n+1)=M(n+1,n+1)+H(2,2);



%partII L2-norm
H1=[h/3, h/6; h/6, h/3];

M1=zeros(n+1,n+1);

for i=1:n
    for j=i:i+1
        for k=i:i+1
            M1(j,k)=M1(j,k)+H1(j-i+1,k-i+1);
        end
    end
end


%% ~~~~~~~~~finished the computation of \|p\|_delta^2~~~~~~~~~~~~~~~~~~~~~

M=M+M1;


%% ~~~~~~~~~ the computation of \|p\|_0^2~~~~~~~~~~~~~~~~~~~~~


%%          
%N for L2-space
N=eye(n,n);
N=h*N;




%%
%SX & SY
SX=chol(M);
%SX=SX';
SY=chol(N);
%SY=SY';


%%
%B
hii=(delta-3*h)/(3*h)*ones(n,1);

H1=diag(hii);
H2=diag(((3*h-2*delta)/(3*h))*ones(n-1,1),-1);
H3=diag((delta/(3*h))*ones(n-2,1),-2);

B=H1+H2+H3;
B=[B;zeros(1,n-2),delta/(3*h),(3*h-2*delta)/(3*h)];
B(1,n)=delta/(3*h);

Sigma=inv(SY)'*B'*inv(SX);
s=svd(Sigma);
s_min=min(s);
W=s_min;

%end


