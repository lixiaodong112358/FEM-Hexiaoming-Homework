%%%运行求解二维非稳态NS方程的求解程序，2021/8/28，李晓东
clc
clear
close all

%%
left=0;
right=1;
bottom=-0.25;
top=0;
N1=6;                       %%%水平方向的网格数目
N2=6;                       %%%垂直方向的网格数目
h1=(right-left)/N1;
h2=(top-bottom)/N2;
%%  时间步长不能太大，不然压力的结果与解析解对不上  2021/9/14
Start=0;
End=0.005;
dt=0.005;
%%
basis_type_u=202;
basis_type_p=201;
% basis_der_x_y_trial=[1 0;0 1];
% basis_der_x_y_test=[1 0;0 1];
basis_der_x_y_test_b=[0 0];
s=[1 1];                    %%误差阶数
%% Solve
[error_u,error_p,solution]=FEM_solver_2D_Unsteady_NS(Start,End,dt,left,right,bottom,top,h1,h2,basis_type_u,basis_type_p,basis_der_x_y_test_b,s);


%% figure
% plot(solution)