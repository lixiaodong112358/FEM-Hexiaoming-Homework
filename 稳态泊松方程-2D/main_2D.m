%%%运行求解二维稳态泊松方程的求解程序，2021/5/14，李晓东
%%%注意这个程序 基本调试通过 如果有什么问题 请检查边界条件
clc
clear
close

left=-1;
right=1;
bottom=-1;
top=1;
N1=10;                       %%%水平方向的网格数目
N2=10;                       %%%垂直方向的网格数目
h1=(right-left)/N1;
h2=(top-bottom)/N2;
basis_type_trial=201;
basis_type_test=201;
basis_der_x_y_trial=[1 0;0 1];
basis_der_x_y_test=[1 0;0 1];
basis_der_x_y_test_b=[0 0];
s=[0 0];
%% Solve
[error,solution]=FEM_solver_2D_Poisson(left,right,bottom,top,h1,h2,basis_type_trial,basis_der_x_y_trial,basis_type_test,basis_der_x_y_test,basis_der_x_y_test_b,s);


%% figure
plot(solution)