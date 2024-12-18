%%%��������ά����̬Stokes���̵�������2021/9/7��������
clc
clear
close all
%% ʱ��
Start=0;
End=0.1;
dt=0.01;
theta=0.5; 
%% �ռ�
left=0;
right=1;
bottom=-0.25;
top=0;                 
N1=8;                       %%%ˮƽ�����������Ŀ
N2=8;                       %%%��ֱ�����������Ŀ
h1=(right-left)/N1;
h2=(top-bottom)/N2;
%%
basis_type_u=202;
basis_type_p=201;
basis_der_x_y_trial=[1 0;0 1];
basis_der_x_y_test=[1 0;0 1];
basis_der_x_y_test_b=[0 0];
s=[0 0];                    %%������
%% Solve
[error_u,error_p,solution]=FEM_solver_2D_Unsteady_stokes(theta,Start,End,dt,left,right,bottom,top,h1,h2,basis_type_u,basis_der_x_y_trial,basis_type_p,basis_der_x_y_test,basis_der_x_y_test_b,s);

error_u

error_p
%% figure
% plot(solution)