%%%��������ά���Է��̵�������2021/7/4��������
clc
clear
close

left=0;
right=1;
bottom=0;
top=1;
N1=10;                       %%%ˮƽ�����������Ŀ
N2=10;                       %%%��ֱ�����������Ŀ
h1=(right-left)/N1;
h2=(top-bottom)/N2;
basis_type_trial=202;
basis_type_test=202;
basis_der_x_y_trial=[1 0;0 1];
basis_der_x_y_test=[1 0;0 1];
basis_der_x_y_test_b=[0 0];
s=[0 0];                    %%������
tic
%% Solve
[error,solution]=FEM_solver_2D_elasticity(left,right,bottom,top,h1,h2,basis_type_trial,basis_der_x_y_trial,basis_type_test,basis_der_x_y_test,basis_der_x_y_test_b,s);
toc

%% figure
% plot(solution)