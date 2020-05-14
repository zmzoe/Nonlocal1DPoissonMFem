

 L=cell(3,2);
 VV=cell(3,1);
 delta=0.5;

for i=1:3
    m=3^(i-1)*5;
    [U,D,V,l]=Bmatrix(delta, m);
    VV{i,1}=V;
    h=delta/m;
    l(size(l))=[];
    L{i,1}=l;
    L{i,2}=l./h*sqrt(3/2);
end

