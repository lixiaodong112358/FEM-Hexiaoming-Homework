function [error,solution]=FEM_solver_1D_Poisson(left,right,h,basis_type_trial,basis_der_x_trial,basis_type_test,basis_der_x_test,basis_der_x_test_b,Gauss_type,s)
% basis_type_trial==101:1D linear
% N:单元个数
% s:误差函数的导数阶
%%%经过测试right值取得很大的话，矩阵求解会出现奇异值
%%%针对一维的泊松方程
%%%经过测试本求解器误差达到了课件的要求；
%%%经过测试，本求解器达到了二次基函数时的误差要求，2021/4/27，李晓东；
%%
% clear
% clc
% close
% left=0;
% right=5;
% h=0.05;  %%%时间步长的选取，暂时不能乱写，因为你必须保证P矩阵的右端点还是1，随意选取的话，MATLAB生成的向量不能保证右端点还是1
% 因此，为了保证不出现此类问题，空间步长暂时取类似于1/12,1/24,1/48这种  李晓东/2021/4/26；
% 这种前处理的工作暂时先放一放;
% basis_type_trial=101;
% basis_type_test=101;
% basis_der_x_trial=1;
% basis_der_x_test=1;
% basis_der_x_test_b=0;
% Gauss_type=4;   %高斯节点取四个
%%

[P,T]=generate_PbTb(left,right,h,101);
[Pb,Tb]=generate_PbTb(left,right,h,basis_type_trial);
N=size(P,2);
Nb=size(Pb,2);
if basis_type_trial==101
  Pb_trial=P;
  Tb_trial=T;
  number_of_local_basis_fun_trial=2;
elseif basis_type_trial==102
  Pb_trial=Pb;
  Tb_trial=Tb;
  number_of_local_basis_fun_trial=3;
end
if basis_type_test==101
  Pb_test=P;
  Tb_test=T;
  number_of_local_basis_fun_test=2;
elseif basis_type_test==102
  Pb_test=Pb;
  Tb_test=Tb;
  number_of_local_basis_fun_test=3;
end
%%%边界节点
boundary_type='DD';
boundarynodes=generate_boundarynodes(Pb,boundary_type);
%% 组装系数矩阵A
matrix_size=[max(size(Pb)) max(size(Pb))];
number_of_elements=size(T,2);
A=assemble_matrix_1D('function_c',Gauss_type,matrix_size,number_of_elements,P,T,Tb_trial,Tb_test,number_of_local_basis_fun_trial,number_of_local_basis_fun_test,basis_type_trial,basis_der_x_trial,basis_type_test,basis_der_x_test);

%%%%%组装向量b
%%
b=assemble_vector_1D('function_f',P,T,Tb,Nb,number_of_elements,number_of_local_basis_fun_test,Gauss_type,basis_type_test,basis_der_x_test_b);
%%%%边界条件的处理
%%
[A0,b]=treat_Dirichlet_boundary('function_g',Pb,A,b,boundarynodes);
%%

% tic
% solution=pinv(full(A))*b;
solution=A0\b;

% solution=linsolve(full(A0),b);
% toc
% plot(P,solution,'*',P,P.*cos(P),'o')
%%%%计算误差
error=compute_Hs_error(P,T,Tb,s,solution,Gauss_type,'fx',number_of_local_basis_fun_test,basis_type_test);
end
% error=compute_max_FE_nodes(?);
% end