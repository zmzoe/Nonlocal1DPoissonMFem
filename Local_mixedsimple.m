

% Mixed finite element method 
% 1D
% P1-P0
% simple version

%-----------------------------
% local
%-----------------------------

% mesh size
N = 4;
h = 1 / N;

% true solution and RHS
u = @(x) sin(2*pi*x);
p = @(x) 2*pi*cos(2*pi*x);
f = @(x) 4*pi^2*sin(2*pi*x);




% FEM matrix
A = spdiags([ones(N,1),4*ones(N,1),ones(N,1)], -1:1, N, N);
A(1,N) = 1;
A(N,1) = 1;
A = A * h/6;

B = spdiags([ones(N,1), -ones(N,1)],-1:0,N,N);
B(1,N)=1;

F = -1*[zeros(N,1);Fgauss(h,f)];

M=[A,B;B',zeros(N,N)];


   
% FEM solution
xx_h = gmres(M, F, [], 1e-12);

p_h = [xx_h(1:N);xx_h(1)];
u_h = xx_h(N+1:2*N);

% interpolated solution
p_r = p(linspace(0, 1, N+1)');
u_r = u(linspace(0, 1-h, N)');


% error
err_p = p_h - p_r;
err_u = u_h - u_r;

max(abs(err_p))
max(abs(err_u))

% figure
figure
    plot(p_h,'LineWidth',2);
    hold on
    plot(p_r);
    legend('numP','realP')
    
    
   figure 
    for i=1:N
        sx=[(i-1)*h,i*h];
        sy=[u_h(i),u_h(i)];
        sw=[u(sx(1)),u(sx(2))]-sy;
        o1=line(sx,sy,'LineWidth',2);
        hold on
        o3=line(sx,sw,'Color',[1 0.8 1],'LineWidth',2);
        hold on
    end
    hold on
    o2=plot(0:h:1,u(0:h:1),'Color',[1 0.5 0],'LineWidth',1);
    legend([o1,o2,o3],{'numerical-u', 'real-u','errU'})
    xlabel('x')
    ylabel('function value')
    title('NumSol-U')
    
    
    
    
   
    
    
    


