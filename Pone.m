%
% Frame
% Import alpha_i, beta_i .   row vector
% Output  (Galpha,beta)

function [GAlphaBeta]=Pone(alpha,beta, delta, m )
%%  Parameter

h = delta/m;
n = 1/h;

%% G alpha

  
B_odd  = sparse(2*n,2*n);
B_even = sparse(2*n,2*n);

vv_odd  = zeros(m+1,1);
vv_even = zeros(m+1,1);

uu_odd  = zeros(m+1,1);
uu_even = zeros(m+1,1);

J_odd  = zeros(m+1,1);
J_even = zeros(m+1,1);

if m == 1
    B(1,2*n-1) = h^2/8;
    B(1,2*n) = 5*h^2/24;
    B(1,1) = -5*h^2/24;
    B(1,2) = -h^2/8;
    
    for i=2:n
        B(2*i-1,2*i-3) = h^2/8;
        B(2*i-1,2*i-2) = 5*h^2/24;
        B(2*i-1,2*i-1) = -5*h^2/24;
        B(2*i-1,2*i) = -h^2/8;
    end
    
    B(2,2*n-1) = h^2/24;
    B(2,2*n) = h^2/8;
    B(2,1) = h^2/24;
    B(2,2) = -5*h^2/24;
    
    for i=2:n
        B(2*i,2*i-3) = h^2/24;
        B(2*i,2*i-2) = h^2/8;
        B(2*i,2*i-1) = h^2/24;
        B(2*i,2*i) = -5*h^2/24;
    end
    
else
   for i = 1:n
    I_odd  = (2*i-1)*ones(m+1,1);
    I_even =  2*i*ones(m+1,1);
    
        if i <= m
        J = [n-m+i : n, 1 : i];
        else
        J = i-m : i;
        end
        
        for j=1 : m+1
        J_odd(j)  = 2*J(j)-1;
        J_even(j) = 2*J(j);
        end
    
       
        vv_odd(1)   =  h^2/8;
        vv_odd(2:m) =  h^2/4;
        vv_odd(m+1) =  (3*h-8*delta)*h/24; 
        vv_even(1)   =  5*h^2/24;
        vv_even(2:m) =  h^2/4;
        vv_even(m+1) =  (h-4*delta)*h/24;
        
        uu_odd(1)   =  h^2/24;
        uu_odd(2:m) =  h^2/4;
        uu_odd(m+1) =  (5*h-4*delta)*h/24;
        uu_even(1)   =  h^2/8;
        uu_even(2:m) =  h^2/4;
        uu_even(m+1) =  (3*h-8*delta)*h/24;
    
       
    B_odd = B_odd + sparse(I_odd ,J_odd ,vv_odd ,2*n,2*n); % odd row odd column 
    B_odd = B_odd + sparse(I_odd ,J_even ,vv_even ,2*n,2*n); % odd row even column 
       
    
    B_even = B_even + sparse(I_even ,J_odd ,uu_odd ,2*n,2*n); % even row odd column  
    B_even = B_even + sparse(I_even ,J_even ,uu_even ,2*n,2*n); % odd row even column 
    
   end  
   B = B_odd + B_even;
end

B = 2/(delta^2)*B;

GAlphaBeta = alpha*B*beta';
end