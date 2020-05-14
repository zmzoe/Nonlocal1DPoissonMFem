

% W is dft matrix W^T=W
% s in (0,delta)=(0,mh)
% h is meshsize
% k is the gs corresponding to the k_th lambda 

% Aa(j,k): eta_k---coefficient of s^2 for piecewise quadratic function over
% I_j



function [Aa,Bb,Cc]=gs4allk(h,m)

n=1/h;
W=dftmtx(n);
W=[W,W(:,1)];
Aa=zeros(m,n);
Bb=zeros(m,n);
Cc=zeros(m,n);




for k=1:n
    for jj=1:m
        j=jj-1;
    Aa(jj,k) = W(k,n-j+1)-2*W(k,n-j)+W(k,n-j-1);
    Bb(jj,k) = 2*(-W(k,n-j+1)*(j+1)+W(k,n-j)*(2*j+1)-W(k,(n-j-1))*j)*h;
    Cc(jj,k) = (W(k,n-j+1)*(j+1)^2+W(k,n-j)*(1-2*j-2*j^2)+W(k,n-j-1)*j^2)*h^2;    
    end
end



end