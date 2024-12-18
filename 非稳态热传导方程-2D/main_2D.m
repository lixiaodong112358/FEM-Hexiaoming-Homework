%%%��������ά��̬���ɷ��̵�������2021/5/14��������
clc
clear
close

left=0;
right=2;
bottom=0;
top=1;
N1=10;                       %%%ˮƽ�����������Ŀ
N2=10;                      %%%��ֱ�����������Ŀ
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
s=[1 1];                    %%������
theta=0.5;
%% Solve
tic
[error,solution]=FEM_solver_2D_heat(theta,Start,End,dt,left,right,bottom,top,h1,h2,basis_type_trial,basis_der_x_y_trial,basis_type_test,basis_der_x_y_test,basis_der_x_y_test_b,s);
toc

%% figure
% plot(solution)