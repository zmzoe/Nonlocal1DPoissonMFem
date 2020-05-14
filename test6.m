


[p,e,t] = initmesh('squareg','hmax',0.1); % mesh 
x = p(1,:); y = p(2,:); % node coordinates 
pif = x.*y; % nodal values of interpolant 
pdesurf(p,t,pif) % plot interpolant
