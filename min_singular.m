
clc


%M=zeros(20,1);

%for i=1:20
    
n=4;

format long

h=1/n;

sy=ones(n,1);
S_y=sqrt(h)*diag(sy);

sx=ones(n+1,1);
S_x=(1/sqrt(h))*(-diag(sx)+diag(sy,1));

B=-diag(sx)+diag(sy,1);
B(n+1,:)=[];
Sigma=inv(S_y)*B*inv(S_x);
s=svd(Sigma);
s_min=min(s);
%M(i,1)=s_min;
%end
%plot(1:20, M)

