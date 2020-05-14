
function res = guasslegendre(fun, a, b, c, d, n)
% -----------------------------------------------------------------------
% Gauss-Legendre��ֵ���ּ�����ػ���
% ����˵��
% fun�����ֱ��ʽ������� fun=@(x,y)f(x,y)
% a,b�����������䣬����
% c,d���ڲ�������䣬������� c=@(x)c(x), d=@(x)d(x)
% n �����ֽ���
% ������
% res�����ֽ��
% -----------------------------------------------------------------------
% 1 ������ֵ�
syms x
p = sym2poly(diff((x^2-1)^(n+1),n+1));
u = roots(p); 

% 2 �������ϵ��
Ak = zeros(n+1,1);
for i=1:n+1
    t = u;
    t(i) = [];
    pn = poly(t);
    fp = @(x)polyval(pn,x)/polyval(pn,u(i));
    Ak(i)=integral(fp,-1,1);
end

% 3 ��������
fa = @(x) (d(x)-c(x))/2.0;
fb = @(x) (d(x)+c(x))/2.0;

% 4 ��е���
[X,Y] = meshgrid(u,u);
[A,B] = meshgrid(Ak,Ak);

U = (b-a)/2.0*X + (b+a)/2.0;
V = fa(U).*Y + fb(U);
w = A.*B.*fa(U).*fun(U,V);

res = (b-a)/2.0*sum(w(:));