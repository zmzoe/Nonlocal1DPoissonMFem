syms t theta S T t

alpha1=@(t) theta*sin(pi/(S*T)*t);
alpha2=@(t) theta*sin(pi/((1-S)*T)*(T-t));


a1=diff(alpha1,t,1);
a2=diff(alpha2,t,1);
