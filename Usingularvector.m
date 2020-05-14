
delta=0.5;
m=5;
h=delta/m;
n=1/h;
f1=@(k,m,n) sqrt(2/n)*cos(2*pi*(k-1)*(m-1)/n);
f2=@(k,m,n) sqrt(2/n)*sin(2*pi*(k-1)*(m-1)/n);

U=zeros(n,1);

U(1,1)=1/sqrt(n);
for i=2:n/2
    U(i,1)=f1(m,i,n);
end

for i=2+n/2:n
    U(i,1)=f2(m,i,n);
end