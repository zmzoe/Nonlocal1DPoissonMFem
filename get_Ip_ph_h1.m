function err=get_Ip_ph_h1(p_h,p_r,h)

n  = 1/h;     

err= 0;
F  = zeros(n,1);

for i=1:n-1
    
    
    F(i)=abs((p_h(i+1)-p_r(i+1))-(p_h(i)-p_r(i)));
  
    err=err+F(i);
end
