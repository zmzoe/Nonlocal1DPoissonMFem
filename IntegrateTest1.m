syms x m delta s h J

f=(x+s-(J-1)*h)/h;
Gf=int(f,s,0,J*h-x)-delta*(x-(J-1)*h)/h;

F=2/(delta^2)*int(Gf,x,(J-1)*h,J*h);               
