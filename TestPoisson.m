%   Main
% 
%   Poisson 不是混合元
%
% 
% continuous p1 p0

%% Parameter

n     =   16;
h     =   1/n;
k     =   1;


%% basis 
%% P1


Alpha = zeros(n, 2*n);

for i = 2:n
    Alpha(i,2*(i-1)) = 1;
    Alpha(i,2*i-1)   = 1;
end

Alpha(1,1)  =  1;
Alpha(1,2*n)=  1;



%% Linear system

A = zeros(n,n);

for i = 1:n
    for j = 1:n
         A(i,j) = Stifnesslocal(Alpha(i,:),Alpha(j,:), n);   
    end
end

Amass = zeros(n,n);

for i = 1:n
    for j = 1:n
         Amass(i,j) = innerproductlocal(Alpha(i,:),Alpha(j,:), n);       
    end
end

%%

ureal= @(x) sin(2*k*pi*x)+cos(2*k*pi*x); 
f=@(x) 4*k^2*pi^2*cos(2*pi*k*x) + 4*k^2*pi^2*sin(2*pi*k*x);

F= FgaussPoissonanybase1(Alpha,f,h);


[xx,flag,relres] = gmres(A,F,[],1e-10);
%du= @(x) 2*k*pi*cos(2*pi*k*x) + 4*pi*cos(2*pi*x).*cos(4*pi*x) - 2*pi*sin(2*pi*x).*sin(4*pi*x);

du=@(x) 2*k*pi*cos(2*pi*k*x) - 2*k*pi*sin(2*pi*k*x);




uh=[xx;xx(1)];%n+1 piecewise linear
err_P_H1=getH1Error_localPoisson_u(uh,du,h);
err_P_H1=sqrt(err_P_H1);

err_P_L2=getL2Error_local_p(uh,ureal,h);
err_P_L2=sqrt(err_P_L2);

%%
err_Iu_u_h_H1=getH1Error_localPoisson_Iu_u_h(uh,ureal,h);
err_Iu_u_h_H1=sqrt(err_Iu_u_h_H1);

q=err_Iu_u_h_H1;




coe=zeros(n,1);
for i=1:n
    coe(i,1)=ureal((i-1)*h)-uh(i);   
end

err_infty=max(coe);

ERR=sqrt(coe'*Amass*coe);

err_P_L2_inter=getL2ErrorInter(uh,ureal,h);
err_P_L2_inter=sqrt(err_P_L2_inter);


err_P_u_Ihu=getL2Error_local_tem(ureal,h);
err_P_u_Ihu=sqrt(err_P_u_Ihu);

  


figure
   plot((0:h:1),uh,'LineWidth',2);
   hold on
   plot((0:h:1),ureal(0:h:1),'--d');
   legend('uh','ureal');
