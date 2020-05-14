function huitu(h,erruIh,errp,errL2u,errH1u,N)
%%%% Plot convergence rates
figure(2);
subplot(1,2,1);
showrateh2(h,erruIh,1,'-*','|u_I-u_h|_1',...
           h,errp,1,'m-+','|| p-p_h||');

fprintf('\n');
disp('Table: Error')
colname = {'#Dof','h','|u_I-u_h|_1','||p-p_h||'};
disptable(colname,N,[],h,'%0.3e',erruIh,'%0.5e',errp,'%0.5e');


figure(2);
subplot(1,2,2);
showrateh3(h,errH1u,2,'-*', '|| Du - Du_h ||', ...
           h,errL2u,2,'k-+', '|| u - u_h ||', ...
           h,erruIh,2,'m-+','|| D u_I - D u_h ||');

fprintf('\n');
disp('Table: Error')
colname = {'#Dof','h','|| u-u_h ||','|| Du-Du_h ||','|| Du_I-Du_h ||'};
disptable(colname,N,[],h,'%0.3e',errL2u,'%0.5e',errH1u,'%0.5e',erruIh,'%0.5e');     
%saveas(gcf,'dane6','png')

end