% mesh size
N = 1000;
h = 1 / N;

% true solution and RHS
u = @(x) sin(2*pi*x);
f = @(x) 4*pi^2*sin(2*pi*x);

% hat function basis
phi = {};
for i = 1 : N
    phi{i} = @(x) 1-abs(x-i*h)/h;
end

% FEM matrix
A = spdiags([-ones(N,1),2*ones(N,1),-ones(N,1)], -1:1, N, N);
A(1,N) = -1;
A(N,1) = -1;
A = A / h;

% FEM RHS
f_h = zeros(N,1);
for i = 1 : N
    f_h(i) = integral(@(x)f(x).*phi{i}(x), (i-1)*h, (i+1)*h, ...
        'AbsTol', 1e-12);
end

% FEM solution
u_h = gmres(A, f_h, [], 1e-12);
% interpolated solution
u_r = u(linspace(h, 1, N)');
% error
err = u_h - u_r;
max(abs(err))

