


delta  =       0.5;
m      =         10;
h      =   delta/m;
n      =       1/h;


%% basis 1 continuous P1- P0
%% P1

Alpha = zeros(1, 2*n);

for i = 2:2:n
    Alpha(1,2*(i-1)) =     -1;
    Alpha(1,2*i-1)   =     -1;
end

for i = 1:2:n-3
    Alpha(1,2*(i+1)) =      1;
    Alpha(1,2*i+3)   =      1;
end

Alpha(1,1)   =      1;
Alpha(1,2*n) =      1;


%% P0

Alpha_p0 = zeros(1, 2*n);

for i = 1:n-1
    Alpha_p0(1,2*i-1)     =  1;
    Alpha_p0(1,2*i)       =  1;
    Alpha_p0(1,2*i+1)     = -1;
    Alpha_p0(1,2*i+2)     = -1;
 
end

%% Linear system


B = zeros(1,1);

for i = 1:1
    for j = 1:1
         B(i,j) = Pone(Alpha(i,:),Alpha_p0(j,:), delta, m );       
    end
end

%aa   =   rank(A);
%aa_1 =   size(A);
%bb   =   rank(B);   
%bb_1 =   size(B);   

