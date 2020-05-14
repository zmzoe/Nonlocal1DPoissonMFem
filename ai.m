 
function [a_i]=ai(X,h,i,x)


a_i=[int(phi_iminus(X,h,i,x)^2,x,X(i),X(i+1)),                   int(phi_iminus(X,h,i,x)*phi_iplus(X,h,i+1,x),x,X(i),X(i+1));
     int(phi_iminus(X,h,i,x)*phi_iplus(X,h,i+1,x),x,X(i),X(i+1)),int(phi_iplus(X,h,i+1,x)^2,x,X(i),X(i+1))];
end