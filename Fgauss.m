

function F=Fgauss(h,f)

XX = 0:h:1;   
n = length(XX)-1;
F=zeros(n,1);
for i=1:n
    F(i) = integral(@(x)f(x), XX(i), XX(i+1), ...
        'AbsTol', 1e-12);
end
