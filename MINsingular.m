%function 
%input n delta
%output W

function [W]=MINsingular(n,delta)
h=1/n;
%delta=1/4*h;%change
%basic information & variables

%%
%M for delta-space
%Part 1 delta-norm-nonlocal
H=zeros(3,3);
H(1,1:3)=[1/5, -1/15, -2/15];
H(2,2:3)=[7/15, -2/5];
H(3,3)=8/15;

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

M(1,1)=M(1,1)+8/15;

M(1,n)=M(1,n)-2/15;
M(n,1)=M(n,1)-2/15;

M(1,n+1)=M(1,n+1)+2/5;
M(n+1,1)=M(n+1,1)+2/5;

M(n,n)=M(n,n)+1/5;

M(n,n+1)=M(n,n+1)-1/15;
M(n+1,n)=M(n+1,n)-1/15;

M(n+1,n+1)=M(n+1,n+1)+8/15;
M=delta/(h^2)*M;



%%
%partII L2-norm
H1=[1, 1/2; 1/2, 1];

M1=zeros(n+1,n+1);

for i=1:n
    for j=i:i+1
        for k=i:i+1
            M1(j,k)=M1(j,k)+H1(j-i+1,k-i+1);
        end
    end
end
M1=h/3*M1;


%%
%partIII delta-norm-noreaction
H2=[1, -1; -1, 1];

M2=zeros(n+1,n+1);

for i=1:n
    for j=i:i+1
        for k=i:i+1
            M2(j,k)=M2(j,k)+H2(j-i+1,k-i+1);
        end
    end
end
M2=(h-delta)/(h^2)*M2;
%
M=M+M1+M2;


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
h1=-1*ones(n+1,1);
h12=ones(n,1);
B1=diag(h1);
B1=B1+diag(h12,-1);
B1(:,n+1)=[];

B1=(h-delta)/h*B1;

B2=diag(h1);
B2=diag(-1*h12, -1)+B2;
B2=B2+diag(2*h12(1:n-1),-2);
B2(:,n+1)=[];
B2(1,n)=2;

B2=delta/(3*h)*B2;

B=B1+B2;

%Sigma=inv(SX)*B*inv(SY);
%Sigma=SX*B*SY;
%Sigma=SY*B'*SX;
%
Sigma=inv(SY)'*B'*inv(SX);
s=svd(Sigma);
s_min=min(s);
W=s_min;

end


