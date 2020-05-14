




delta=0.5;
Sigma=cell(3,1);
for l=1:3
m=5*(3^(l-1));

h=delta/m;
n=1/h;


[A,~,~,~,~,~,~,B]=Bmatrix(delta, m);
Sigma{l,1} = gsvd(B,A);

end