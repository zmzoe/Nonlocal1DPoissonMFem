



%g(s)=int_o^1(eta_k(x-s) phi_1(x))

%lambda_k=int_0^delta rho_delta(s)(g(0)-g(s))
%lambda_k is the eigenvalue of B
% B is a circulant matrix and B(ij)=\int_Ij G_delta^* phi_i
%%

syms c1 w h n 



pj=2;
k=2;

s=0;
qua=(w^((n-pj)*(k-1))-2*w^((n-pj-1)*(k-1))+w^((n-pj-2)*(k-1)))*s^2;

lin=2*(-w^((n-pj)*(k-1))*(pj+1)+(2*pj+1)*w^((n-pj-1)*(k-1))-w^((n-pj-2)*(k-1))*pj)*h*s;

const=(w^((k-1)*(n-pj))*(pj+1)^2+w^((n-pj-1)*(k-1))*(1-2*pj-2*pj^2)+w^((n-pj-2)*(k-1))*pj^2)*h^2;

gs=1/(2*h)*(qua+lin+const);

