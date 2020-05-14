%Test
function [A,L,U,D,V,l1,l,h,B]=Bmatrix(delta, m)
format long


%% Parameter
h      =   delta/m;
n      =   1/h;




%% bases con
Alpha = zeros(n, 2*n);

for i = 2:n
    Alpha(i,2*(i-1)) = 1;
    Alpha(i,2*i-1)   = 1;
end

Alpha(1,1)  =  1;
Alpha(1,2*n)=  1;

A = zeros(n,n);

for i = 1:n
    for j = 1:n
         A(i,j) = Stifnesslocal(Alpha(i,:),Alpha(j,:), n);   
    end
end

Alpha_p0 = zeros(n, 2*n);

for i = 1:n
    Alpha_p0(i,2*i-1)   = 1;
    Alpha_p0(i,2*i)     = 1;
 
end

B = zeros(n,n);

for i = 1:n
    for j = 1:n
         B(i,j) = Pone(Alpha(i,:),Alpha_p0(j,:), delta, m );       
    end
end

%A=B'*B;
%lambda=eig(A);

[U,D,V]=svd(B);%B=UDV'
l=svd(B);
l1=eig(B);
l1=abs(l1);
L=l/h;
end

