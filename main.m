%main
clear;
clc;
format long

%n
W=zeros(100,3);

 

for j=1:3
    delta=10^(-3-j+1);
     for i=3:100
        
         n=i;
       [pS,Sigma]=Mtest4(delta, n);
        %[Sigma]=P1singular(n,delta);
        W(i,j)=Sigma;
    end

end

figure
plot(3:100, W(3:100, 1),'--s','LineWidth',1);
hold on 
plot(3:100, W(3:100, 2),'LineWidth',2);
hold on 
plot(3:8:100, W(3:8:100, 3),'--','LineWidth',2);

set(gca, 'LineWidth', 1.5);
title('Relationship of MinSingularvalues and h')
xlabel('#points(1/h)')
ylabel('MinSigularvalue')
legend('/delta=10^{-3}','/delta=10^{-4}','/delta=10^{-5}')

