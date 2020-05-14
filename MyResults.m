
load dispdata
%results---fix m
par_h=zeros(4,2);
M_re=zeros(4,6); %第k行代表第k层上的各种误差
for k=1:4
    m=5;
    delta=0.4/(2^k);
    par_h(k,1)=delta./m;
    [err_P_delta,err_P_L2,err_U,ERR_u0,ERR0,ERR1]=myNonlocalMixed1(delta,m);
    M_re(k,1)=err_P_delta;
    M_re(k,2)=err_P_L2;
    M_re(k,3)=err_U;
    M_re(k,4)=ERR0;
    M_re(k,5)=ERR1;
    M_re(k,6)=ERR_u0;
end

%results---fix m
D_re=zeros(4,6);
for k=1:4
    delta=0.5;
    m=5*2^k;
    par_h(k,2)=delta./m;
    [err_P_delta,err_P_L2,err_U,ERR_u0,ERR0,ERR1]=myNonlocalMixed1(delta,m);
    D_re(k,1)=err_P_delta;
    D_re(k,2)=err_P_L2;
    D_re(k,3)=err_U;
    D_re(k,4)=ERR0;
    D_re(k,5)=ERR1;
    D_re(k,6)=ERR_u0;
end

Rate_m=zeros(3,6);
Rate_d=zeros(3,6);

for k=1:3
    Rate_m(k,:)=M_re(k,:)./M_re(k+1,:);
    Rate_d(k,:)=D_re(k,:)./D_re(k+1,:);
end

Rate_mm=zeros(4,6);
Rate_dd=zeros(4,6);
Rate_mm(2:4,:)=Rate_m;
Rate_dd(2:4,:)=Rate_d;

ConvergenceRateFigure(D_re,M_re,par_h,Rate_m,Rate_d,Rate_mm,Rate_dd);



function ConvergenceRateFigure(D_re,M_re,par_h,Rate_m,Rate_d,Rate_mm,Rate_dd)

figure(1);
subplot(2,2,1);
p=plot(1:3,log2(Rate_m(:,1)),'-*',1:3,log2(Rate_m(:,2)),'m-+',1:3,log2(Rate_m(:,3)),'k-s');
hold on   
p(1).LineWidth=2;
p(2).LineWidth=2;
p(3).LineWidth=2;
axis([1 3 0 3])
legend('|p-p_h|_{delta}','|| p-p_h||','|| u-u_h||') 
title('optimal-Fixedm')

fprintf('\n');
disp('Table: Error4Fixedm')
colname = {'h','|p-p_h|_{delta}','order','|| p-p_h||','order','|| u-u_h||','order'};
disptable(colname,par_h(:,1),[],M_re(:,1),'%0.5e',Rate_mm(:,1),'%0.3g',M_re(:,2),'%0.5e',Rate_mm(:,2),'%0.3g',M_re(:,3),'%0.5e',Rate_mm(:,3),'%0.3g');

figure(1);
subplot(2,2,2);
p=plot(1:3,log2(Rate_d(:,1)),'-*',1:3,log2(Rate_d(:,2)),'m-+',1:3,log2(Rate_d(:,3)),'k-s');
hold on   
p(1).LineWidth=2;
p(2).LineWidth=1;
p(3).LineWidth=2;
axis([1 3 0 3])
legend('|p-p_h|_{delta}','|| p-p_h||','|| u-u_h||') 
title('optimal-Fixeddelta')


fprintf('\n');
disp('Table: Error4Fixeddelta')
colname = {'h','|p-p_h|_{delta}','order','|| p-p_h||','order','|| u-u_h||','order'};
disptable(colname,par_h(:,2),[],D_re(:,1),'%0.5e',Rate_dd(:,1),'%0.3g',D_re(:,2),'%0.5e',Rate_dd(:,2),'%0.3g',D_re(:,3),'%0.5e',Rate_dd(:,3),'%0.3g');
 

figure(1);
subplot(2,2,4);
p=plot(1:3,log2(Rate_d(:,5)),'-*',1:3,log2(Rate_d(:,4)),'m-+');
hold on   
p(1).LineWidth=2;
p(2).LineWidth=1;
axis([1 3 0 3])
legend('|I_hp-p_h|_{delta}','|| I_hp-p_h||') 
title('superconvergent-Fixeddelta')

figure(1);
subplot(2,2,3);
p=plot(1:3,log2(Rate_m(:,5)),'-*',1:3,log2(Rate_m(:,4)),'m-+');
hold on   
p(1).LineWidth=2;
p(2).LineWidth=2;
axis([1 3 0 3])
legend('|I_hp-p_h|_{delta}','|| I_hp-p_h||') 
title('superconvergent-Fixedm')

fprintf('\n');
disp('Table: Error4Fixedm')
colname = {'h','|I_hp-p_h|_{delta}','|| I_hp-p_h||','|| I_0u-u_h||'};
disptable(colname,par_h(:,1),[],M_re(:,4),'%0.5e',Rate_mm(:,4),'%0.3g',M_re(:,5),'%0.5e',Rate_mm(:,5),'%0.3g',M_re(:,6),'%0.5e',Rate_mm(:,6),'%0.3g');
fprintf('\n');
disp('Table: Error4Fixeddelta')
colname = {'h','|I_hp-p_h|_{delta}','|| I_hp-p_h||','|| I0_h u-u_h||'};
disptable(colname,par_h(:,2),[],D_re(:,4),'%0.5e',Rate_dd(:,4),'%0.3g',D_re(:,5),'%0.5e',Rate_dd(:,5),'%0.3g',D_re(:,6),'%0.5e',Rate_dd(:,6),'%0.3g');
end
