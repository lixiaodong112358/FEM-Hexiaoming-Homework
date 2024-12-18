clear
clc
close
left=0;
right=1.0;
h=1/10;
xx=(left:h:right)';
basis_type_trial=102;
basis_type_test=102;
basis_der_x_trial=1;
basis_der_x_test=1;
basis_der_x_test_b=0;
s=0;            %误差函数导数阶取值0
Gauss_type=4;   %高斯节点取四个
tic
[error,solution]=FEM_solver_1D_Poisson(left,right,h,basis_type_trial,basis_der_x_trial,basis_type_test,basis_der_x_test,basis_der_x_test_b,Gauss_type,s);
toc
% error=max(abs(solution-xx.*cos(xx)))
%% figure
plot(xx,xx.*cos(xx),'--',linspace(left,right,length(solution)),solution,'*');