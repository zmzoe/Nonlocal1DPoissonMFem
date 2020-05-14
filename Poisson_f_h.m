% mesh size
N = 16;
h = 1 / N;
k = 1;

% true solution and RHS
u = @(x) sin(2*k*pi*x).*cos(2*k*pi*x);
f = @(x) 16*k^2*pi^2*sin(2*pi*k*x).*cos(2*pi*k*x);

% hat function basis
phi = {};
for i = 1 : N+1
    phi{i} = @(x) 1-abs(x-(i-1)*h)/h;
end

% FEM matrix
A = spdiags([-ones(N,1),2*ones(N,1),-ones(N,1)], -1:1, N, N);
A(1,N) = -1;
A(N,1) = -1;
A = A / h;

% FEM RHS
f_h = zeros(N,1);
for i = 2 : N
    f_h(i) = integral(@(x)f(x).*phi{i}(x), (i-2)*h, i*h, ...
        'AbsTol', 1e-12);
end

f_h(1)=integral(@(x)f(x).*phi{1}(x),0,h,'AbsTol', 1e-12)+integral(@(x)f(x).*phi{N+1}(x),(N-1)*h,N*h,'AbsTol', 1e-12);

% FEM solution
u_h = gmres(A, f_h, [], 1e-12);
% interpolated solution
u_r = u(linspace(0, 1-h, N)');
% error
err = u_h - u_r;
max(abs(err))

