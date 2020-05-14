

delta  =  0.125;
m  =  5;
h  =  delta/m;
n  =  1/h;



    


%  ------------------------------------------------------------------------
%  remember to add at the end of the computation 2/delta^2 
%  ------------------------------------------------------------------------

M_q=zeros(n,n);     % quadratic term
M_l=zeros(n,n);     % linear term
M_c=zeros(n,n);     % constant term

vv_q=zeros(1,m+2);
vv_l=zeros(1,m+2);
vv_c=zeros(1,m+2);

for i=1:m+1
    II=i*ones(1,m+2);
    JJ=[n-m+i-1:n,1:i];
    
    vv_q(1)=1/(2*h);
    vv_l(1)=m-n+2-i;
    vv_c(1)=h/2*(m-n+2-i);
    
    vv_q(2)=-1/(2*h);
    vv_l(2)=m-n+i;
    vv_c(2)=h*((m-n)*(2*i-m+n)+2-i^2)/2;
    
    for j=3:m
        vv_c(j)=h;
    end   
    
  
    if i==1
    vv_q(m+1)=-1/(2*h);
    vv_l(m+1)=n-m-1;
    vv_c(m+1)=(n+m*(n-1)-n^2/2+1/2)*h;
    else
        vv_q(m+1)=-1/(2*h);
        vv_l(m+1)=i-m-2;
        vv_c(m+1)=h*(2*m*(i-2)-(i-1)*(i-3)+1)/2;
    end
    
    
    vv_q(m+2)=1/(2*h);
    vv_l(m+2)=m-i;
    vv_c(m+2)=h*(i^2/2 - m*i);
    
    M_q = M_q + sparse(II ,JJ ,vv_q ,n,n);
    M_l = M_l + sparse(II ,JJ ,vv_l ,n,n);
    M_c = M_c + sparse(II ,JJ ,vv_c ,n,n);
end


for i=m+2:n
    II=i*ones(1,m+2);
    JJ=(i-m-1:i);
    
    vv_q(1)=1/(2*h);
    vv_l(1)=m+2-i;
    vv_c(1)=h/2*(m+2-i)^2;
    
    vv_q(2)=-1/(2*h);
    vv_l(2)=i-m;
    vv_c(2)=h*(2-(m-i)^2)/2;
   
    for j=3:m
        vv_c(j)=h;
    end
    
    vv_q(m+1)=-1/(2*h);
    vv_l(m+1)=i-m-2;
    vv_c(m+1)=h*(2*m*(i-2)-(i-1)*(i-3)+1)/2;
    
    vv_q(m+2)=1/(2*h);
    vv_l(m+2)=m-i;
    vv_c(m+2)=h*(i^2/2 - m*i);
    
    M_q = M_q + sparse(II ,JJ ,vv_q ,n,n);
    M_l = M_l + sparse(II ,JJ ,vv_l ,n,n);
    M_c = M_c + sparse(II ,JJ ,vv_c ,n,n);
    
end



    M_q = 2/delta^2*M_q;
    M_l = 2/delta^2*M_l;
    M_c = 2/delta^2*M_c;


    
    
   
%  ========================================================================
%                    G_delta part
%  ========================================================================   
    
    
    
    T=zeros(n,n);
    
    Quartic=zeros(n,n);
    Cubic=zeros(n,n);
    Quadratic=zeros(n,n);
    Linear=zeros(n,n);
    Constant=zeros(n,n);
    
    for i=1:n
        for j=1:n
            for k=1:n
            Quartic(i,j)=M_q(i,k)*M_q(j,k);
            Cubic(i,j)=M_q(i,k)*M_l(j,k)+M_q(j,k)*M_l(i,k);
            Quadratic(i,j)=M_q(i,k)*M_c(j,k)+M_q(j,k)*M_c(i,k)+M_l(i,k)*M_l(j,k);
            Linear(i,j)=M_l(i,k)*M_c(j,k)+M_l(j,k)*M_c(i,k);
            Constant(i,j)=M_c(i,k)*M_c(j,k);
            
            T(i,j)=T(i,j)+Quartic(i,j)*(k^5/5 - (k - 1)^5/5)*h^5;
            T(i,j)=T(i,j)+Cubic(i,j)*(k^4/4 - (k - 1)^4/4)*h^4;
            T(i,j)=T(i,j)+Quadratic(i,j)*(k^2 - k + 1/3)*h^3;
            T(i,j)=T(i,j)+Linear(i,j)*(k - 1/2)*h^2;
            T(i,j)=T(i,j)+Constant(i,j)*h;
            end
        end
    end
    
    
    
%  ========================================================================
%                    L2 part
%  ========================================================================
    

    T0=zeros(n,n);
    
    for i=2:n-1
        T0(i,i)=2*h/3;
        T0(i,i-1)=h/6;
        T0(i-1,i)=T0(i,i-1);
    end
    
    
    T0(1,1)=2*h/3;
    T0(n,n)=2*h/3;
    T0(1,n)=h/6;
    T0(n,1)=h/6;
    T0(n-1,n)=h/6;
    T0(n,n-1)=h/6;
    
    TT=T0+T;
    
    
    
 lambda=eig(TT);
    