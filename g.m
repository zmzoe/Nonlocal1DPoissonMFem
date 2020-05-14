function g_i=g(P,i,X,h,x)

k=size(P,2);%n+1=k

[phi_m]=phi_iminus(X,h,i,x);
[phi_p]=phi_iplus(X,h,i+1,x);

if i<k
g_i=P(i)*phi_m+P(i+1)*phi_p;
else
g_i=P(i-k+1)*phi_m+P(i-k+2)*phi_p;  
end

end