

format short

%% MAIN


delta=0.5;
syms x s 


    m=5;
    [D]=SingularValue(delta,m,x,s);
    







%{
figure
plot([2^3,2^4,2^5,2^6], WWW(1,:),'LineWidth',2);
hold on
plot([2^3,2^4,2^5,2^6], WWW(2,:),'LineWidth',2);
hold on
plot([2^3,2^4,2^5,2^6], WWW(3,:),'LineWidth',2);
hold on
plot([2^3,2^4,2^5,2^6], WWW(4,:),'LineWidth',2);
hold on
plot([2^3,2^4,2^5,2^6], WWW(5,:),'LineWidth',2);


legend('MSin','MEigXX','MinSinB','MinEigY','MinEigDelta')
set(gca, 'LineWidth', 1.5);
title('Relationship of MinEigvalues and h')
xlabel('#points(1/h)')
ylabel('MinEigvalue')



figure
plot([2^3,2^4,2^5,2^6], WWW(2,:),'--','LineWidth',2);

set(gca, 'LineWidth', 1.5);
title('Relationship of MinEigvalues and h')
xlabel('#points(1/h)')
ylabel('MinEigvalue')

figure
plot([2^3,2^4,2^5,2^6], WWW(1,:),'--','LineWidth',2);
hold on
plot([2^3,2^4,2^5,2^6], WWW(2,:),'LineWidth',2);
legend('MinSin','MinEigOFD')
set(gca, 'LineWidth', 1.5);
title('Relationship of MinEigvalues and h')
xlabel('#points(1/h)')
ylabel('MinEigvalue')

%legend('h=10^{-1}','h=10^{-2}','h=10^{-3}')

%}