


function [ureal,preal]=myrealu(x,f1,f2)

ureal=subs(f1(x),x,x+1).*(x>=-1 & x<-1/2)+f2(x).*(x>=-1/2 & x<0)+f1(x).*(x>=0 & x<1/2)+subs(f2(x),x,x-1).*(x>=1/2 & x<1)+subs(f1(x),x,x-1).*(x>=1 & x<3/2)+subs(f2(x),x,x-2).*(x>=3/2 & x<=2);
preal=diff(ureal,x,1);
end