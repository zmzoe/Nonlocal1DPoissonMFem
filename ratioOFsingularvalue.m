


%l1 l2 l3
m=5;
k=3;
s=1;
delta=0.5;
[U1,~,V1,l1]=Bmatrix(delta, m);
[U2,~,V2,l2]=Bmatrix(delta/s, k*m);
[U3,~,V3,l3]=Bmatrix(delta/(s^2), k^2*m);



n=length(l1);



if mod(m,2)
l1=l1(1:n-1);
l2=l2(1:k*n-1);
l3=l3(1:k^2*n-1);
else
    l1=l1(1:n-2);
    l2=l2(1:k*n-2);
    l3=l3(1:k^2*n-2);
end
   
t=length(l1);
l1=flip(l1);
l2=flip(l2);
l3=flip(l3);

rate12=zeros(k*t,1);

for i=1:t
    for j=1:k
    rate12((i-1)*k+j,1)=l1(i)/l2((i-1)*k+j);
    end
end

rate23=zeros(k^2*t,1);

for i=1:k*t
    for j=1:k
    rate23((i-1)*k+j,1)=l2(i)/l3((i-1)*k+j);
    end
end


