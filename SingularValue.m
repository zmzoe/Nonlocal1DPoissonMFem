%function [D]=SingularValue(delta,m,x, s)

%h=delta/M
%Omega=(0,1) periodical boundary condition

delta=0.5;
syms x s 


    m=5;
h=delta/m;
n=1/h;
X=0:h:(n+m)*h;


%% ~~~~~~~~~~xi < x < xi+1~~~~~~~~~~~~~~~G_delta^* p  

AA=sym(zeros(n,n)); %G_delta(x_i, x_i+1)
D=sym(zeros(n,n)); %int ((G)^2)
%F=sym(sparse(n,n)); %int(G_delta)


tic
for i=1:(n-m-1)
    
    AA(i,i)=int(phi_iminus(X,h,i,x+s),s,0,X(i+1)-x)-int(phi_iminus(X,h,i,x),s,0,delta);
    AA(i,i+1)=int(phi_iplus(X,h,i+1,x+s),s,0,X(i+1)-x)+int(phi_iminus(X,h,i+1,x+s),s,X(i+1)-x,X(i+2)-x)-int(phi_iplus(X,h,i+1,x),s,0,delta);
    AA(i,i+m)=int(phi_iplus(X,h,i+m,x+s),s,X(i+m-1)-x,X(i+m)-x)+int(phi_iminus(X,h,i+m,x+s),s,X(i+m)-x,delta);  
    AA(i,i+m+1)=int(phi_iplus(X,h,i+m+1,x+s),s,X(i+m)-x,delta);
     
    for j=(i+2):(i+m-1)
        
      % A(i,j)=int(phi_iplus(X,h,j,x+s),s,X(j-1)-x,X(j)-x)+int(phi_iminus(X,h,j,x+s),s,X(j)-x,X(j+1)-x);
      AA(i,j)=h;
          
    end
end

for i=(n-m):n
    AA(i,i)=int(phi_iminus(X,h,i,x+s),s,0,X(i+1)-x)-int(phi_iminus(X,h,i,x),s,0,delta);
    
    if i<n
        AA(i,i+1)=int(phi_iplus(X,h,i+1,x+s),s,0,X(i+1)-x)+int(phi_iminus(X,h,i+1,x+s),s,X(i+1)-x,X(i+2)-x)-int(phi_iplus(X,h,i+1,x),s,0,delta);
    else
        AA(i,i+1-n)=AA(i,i+1-n)+int(phi_iplus(X,h,i+1,x+s),s,0,X(i+1)-x)+int(phi_iminus(X,h,i+1,x+s),s,X(i+1)-x,X(i+2)-x)-int(phi_iplus(X,h,i+1,x),s,0,delta);
    end
    
    
    if i<(n-m+1)
        AA(i,i+m)=int(phi_iplus(X,h,i+m,x+s),s,X(i+m-1)-x,X(i+m)-x)+int(phi_iminus(X,h,i+m,x+s),s,X(i+m)-x,delta);
    else
        %AA(i,i+m-n)=AA(i,i+m-n)+int(phi_iplus(X,h,i+m,x+s),s,X(i+m-1)-x,X(i+m)-x)+int(phi_iminus(X,h,i+m,x+s),s,X(i+m)-x,delta);
        AA(i,i+m-n)=int(phi_iplus(X,h,i+m,x+s),s,X(i+m-1)-x,X(i+m)-x)+int(phi_iminus(X,h,i+m,x+s),s,X(i+m)-x,delta);
    end
    
    %AA(i,i+m+1-n)=AA(i,i+m+1-n)+int(phi_iplus(X,h,i+m+1,x+s),s,X(i+m)-x,delta);
    AA(i,i+m+1-n)=int(phi_iplus(X,h,i+m+1,x+s),s,X(i+m)-x,delta);
    
    for j=(i+2):(i+m-1)
        
        %if j< n
         %    A(i,j)=int(phi_iplus(X,h,j,x+s),s,X(j-1)-x,X(j)-x)+int(phi_iminus(X,h,j,x+s),s,X(j)-x,X(j+1)-x);
         %else
         %    A(i,j-n+1)=A(i,j-n+1)+int(phi_iplus(X,h,j,x+s),s,X(j-1)-x,X(j)-x)+int(phi_iminus(X,h,j,x+s),s,X(j)-x,X(j+1)-x);    
        %end   
        if j<=n
            AA(i,j)=h;
         else
             AA(i,j-n)=AA(i,j-n)+h;    
        end   
          
    end
end



AA=2/(delta^2)*AA;
    
  %% ~~~~~~~~~~~\|p\|^2
  
  %{
  AA=sym(zeros(n,n));
  tic
  for i=1:(n-1)
      a_i=ai(X,h,i,x);
      
      AA(i,i)=A(i,i)+a_i(1,1);
      AA(i,i+1)=A(i,i+1)+a_i(1,2);
      AA(i+1,i)=A(i+1,i)+a_i(2,1);
      AA(i+1,i+1)=A(i+1,i+1)+a_i(2,2);
     
  end
  toc
  %}

  AAL=zeros(n,n);

  for i=1:(n-1)
      a_i=ai(X,h,i,x);
      for j=i:i+1
          for k=i:i+1
               AAL(j,k)=AAL(j,k)+a_i(j-i+1,k-i+1);
          end
      end
  end


  
  a_n=ai(X,h,n,x);  
  
  AAL(n,n)=AAL(n,n)+a_n(1,1);
  AAL(n,1)=AAL(n,1)+a_n(2,1);
  AAL(1,n)=AAL(1,n)+a_n(1,2);
  AAL(1,1)=AAL(1,1)+a_n(2,2);
  
 
  
                      
              
  %% ~~~~~~~~~~~~~~~~~~~~~~~~\|G_delta p\|^2~~~~~~~~~~~~~~~~~~~~~~



   for i=1:n
        M=AA(i,:)'*AA(i,:);
        M=int(M,x,X(i),X(i+1));
        D=D+M;
        
        
  end
                
                
   toc   
    
%% ~~~~~~~~~~~`\|p\|_delta^2~~~~~~~~~~~~~~~~~
D=double(D);

%ee=eig(D);
%D=D+AA;





%% ~~~~~~~~~~~~~~~~~~~~ F ~~~~~~~~~~~~~~~~~~

FF=sym(zeros(n,n));
F=sym(zeros(n,n));
for i=1:n
    %for j=1:n
    %    FF(i,j)=int(AA(i,j),x,X(i),X(i+1));      
    %end
    FF(i,:)=int(AA(i,:),x,X(i),X(i+1)); 
    F(i,:)=F(i,:)+FF(i,:);
end

F=double(F);



%%

%% ~~~~~~~~~~~~~~~ \|v\|_0^2 ~~~~~~~~~~~~~~~~~~   
     
N=h*eye(n);
     
%% ~~~~~~~~~~~~~~~~~SX & SY ~~~~~~~~~~~~~~~~~~~

%{
tic
SX=chol(D);
SY=chol(N);
toc
%}

%{




eee=eig(N);




TT=inv(SY)'*F*inv(SX);

tic
W = svd(TT);
toc
W(W<=(10^-14))=[];

s_min=min(W);

%}




