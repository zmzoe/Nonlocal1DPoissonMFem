


m=5;
Sigma=cell(3,1);
for l=1:3
delta=0.5/(2^(l-1));

h=delta/m;
n=1/h;
T = diag(2*ones(n,1));
T = T+diag(-ones(n-1,1), -1);
T = T+diag(-ones(n-1,1), 1);
T(1,n) = -1;
T(n,1) = -1;
T = T / h;

[~,~,~,~,~,~,~,B]=Bmatrix(delta, m);
Sigma{l,1} = gsvd(B',T);

end

l1=Sigma{1,1};
l2=Sigma{2,1};
l3=Sigma{3,1};