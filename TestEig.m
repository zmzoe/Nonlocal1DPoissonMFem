

%{

[~,~,~,~,~,~,~,B]=Bmatrix(0.5, 5);
B1=B';
[~,~,~,~,~,~,~,B]=Bmatrix(0.25, 5);
B2=B';


lambda1=eig(B1);
lambda2=eig(B2);
l1=abs(lambda1);
l2=abs(lambda2);
sort(l1);
l1=sort(l1);
l2=sort(l2);


m=5;
delta=0.25;
h=delta/m;
n=1/h;

beta=2*ones(1,m-2);
beta=[1/3,5/3,beta,(5-3*m)/3,(1-3*m)/3];
beta=h^2/(delta^2)*beta;
fullbeta=[beta,zeros(1,n-m-2)];

betahat=fft(fullbeta)';

SigmaB=abs(betahat);

W=dftmtx(n);
betatilde=W*fullbeta';
%}


%syms x m 
%{
m=5;
n=10;
norm=zeros(n,1);
for i=1:n
    x=(i-1)*2*pi/n;
re=2+cos(x)-cos(m*x)*cos(x)-2*cos(m*x)+3*m*sin(m*x)*sin(x);
im=3*m*sin(x)*cos(m*x)+sin(m*x)*cos(x)+2*sin(m*x);

norm(i)=(re^2+im^2);
end


%}


m=5;
delta=0.125/32;
h=delta/m;
n=1/h;

w=@(k) exp(-1i*(2*pi)*(k-1)/n); 


alpha1=@(k) (1-w(k)^m)*(1+4*w(k)+w(k)^2)/(3*(1-w(k)))-m*w(k)^m*(1+w(k));

lambdak=@(k) 1/(m^2)*(alpha1(k));

Lambda=zeros(n,1);

for i=1:n
    Lambda(i,1)=abs(lambdak(i));
end
betaT=n*[-1,2,-1,zeros(1,n-3)];
betaS=1/n*[1,zeros(1,n-1)];
hatT=abs(fft(betaT))';
hatS=abs(fft(betaS))';

BB=diag(Lambda);
TT=diag(hatT);
SS=diag(hatS);
mu=BB./sqrt(TT*SS);
u=diag(mu);


%{
beta=2*ones(1,m-2);
beta=[1/3,5/3,beta,(5-3*m)/3,(1-3*m)/3];
beta=h^2/(delta^2)*beta;
fullbeta=[beta,zeros(1,n-m-2)];

betahat=fft(fullbeta)';

SigmaB=abs(betahat);

W=dftmtx(n);
betatilde=W*fullbeta';
%}

%{
%% n=20
n=40;
betaB=[1,-1,zeros(1,n-2)];
betaT=n*[-1,2,-1,zeros(1,n-3)];
betaS=1/n*[1,zeros(1,n-1)];

hatB=abs(fft(betaB))';

hatT=abs(fft(betaT))';

hatS=abs(fft(betaS))';


BB=diag(hatB);
TT=diag(hatT);
SS=diag(hatS);

mu=BB./sqrt(TT*SS);

u=diag(mu);

%}
