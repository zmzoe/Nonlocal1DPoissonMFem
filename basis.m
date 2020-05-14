function [phi,cons]=basis(X,h,i,x)
if i==1
    phi=@(x) (X(i+1)-x)/h;
else
    if (x>=X(i-1) && x<=X(i))
    phi=@(x) (x-X(i-1))/h;
    else
    if (x>=X(i) && x<=X(i+1))
        phi=@(x) (X(i+1)-x)/h;
    else
        phi=@(x)0;
    end
    end
end

if (x>=X(i) && x<=X(i+1))
    cons=@(x) 1;
else
    cons=@(x) 0;
end

end
