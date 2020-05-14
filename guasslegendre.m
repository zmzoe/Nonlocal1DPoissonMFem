
function res = guasslegendre(fun, a, b, c, d, n)
% -----------------------------------------------------------------------
% Gauss-Legendre数值积分计算二重积分
% 参数说明
% fun：积分表达式函数句柄 fun=@(x,y)f(x,y)
% a,b：外层积分区间，常数
% c,d：内层积分区间，函数句柄 c=@(x)c(x), d=@(x)d(x)
% n ：积分阶数
% 输出结果
% res：积分结果
% -----------------------------------------------------------------------
% 1 计算积分点
syms x
p = sym2poly(diff((x^2-1)^(n+1),n+1));
u = roots(p); 

% 2 计算求积系数
Ak = zeros(n+1,1);
for i=1:n+1
    t = u;
    t(i) = [];
    pn = poly(t);
    fp = @(x)polyval(pn,x)/polyval(pn,u(i));
    Ak(i)=integral(fp,-1,1);
end

% 3 变量代换
fa = @(x) (d(x)-c(x))/2.0;
fb = @(x) (d(x)+c(x))/2.0;

% 4 机械求积
[X,Y] = meshgrid(u,u);
[A,B] = meshgrid(Ak,Ak);

U = (b-a)/2.0*X + (b+a)/2.0;
V = fa(U).*Y + fb(U);
w = A.*B.*fa(U).*fun(U,V);

res = (b-a)/2.0*sum(w(:));