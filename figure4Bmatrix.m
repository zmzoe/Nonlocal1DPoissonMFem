


m=5;
k=2;
delta=0.5;

[U1,~,V1,~]=Bmatrix(delta, m);
[U2,~,V2,~]=Bmatrix(delta/k, m);
[U3,~,V3,~]=Bmatrix(delta/(k^2), m);

[n,~]=size(U1);

figure 
plot(U1(n-1,:))
hold on
plot(U1(n-2,:))

figure 
plot(U2(n-1,:))
hold on
plot(U2(n-2,:))

figure

plot(U3(n-1,:))
hold on
plot(U3(n-2,:))