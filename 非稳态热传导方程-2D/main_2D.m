%%%运行求解二维稳态泊松方程的求解程序，2021/5/14，李晓东
clc
clear
close

left=0;
right=2;
bottom=0;
top=1;
N1=10;                       %%%水平方向的网格数目
N2=10;                      %%%垂直方向的网格数目
h1=(right-left)/N1;
h2=(top-bottom)/N2;
dt=0.1;
Start=0;
End=Start+2*dt;
basis_type_trial=201;
basis_type_test=201;
basis_der_x_y_trial=[1 0;0 1];
basis_der_x_y_test=[1 0;0 1];
basis_der_x_y_test_b=[0 0];
s=[1 1];                    %%误差阶数
theta=0.5;
%% Solve
tic
[error,solution]=FEM_solver_2D_heat(theta,Start,End,dt,left,right,bottom,top,h1,h2,basis_type_trial,basis_der_x_y_trial,basis_type_test,basis_der_x_y_test,basis_der_x_y_test_b,s);
toc

%% figure
% plot(solution)