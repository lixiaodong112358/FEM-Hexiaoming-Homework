%%%��������ά����̬���Է��̵�������2021/9/9��������
clc
clear
close
%%
Start=0;
End=1;
dt=0.1;
%%
left=-1;
right=1;
bottom=-1;
top=1;
N1=4;                       %%%ˮƽ�����������Ŀ
N2=4;                       %%%��ֱ�����������Ŀ
h1=(right-left)/N1;
h2=(top-bottom)/N2;
basis_type_trial=201;
basis_type_test=201;
basis_der_x_y_test_b=[0 0];
s=[0 0];                    %%������
tic
%% Solve
[error,solution]=FEM_solver_2D_Unstable_elasticity(Start,End,dt,left,right,bottom,top,h1,h2,basis_type_trial,basis_type_test,basis_der_x_y_test_b,s);
toc

%% figure
% plot(solution)