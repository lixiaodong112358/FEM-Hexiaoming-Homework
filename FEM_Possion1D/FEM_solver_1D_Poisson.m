function [error,solution]=FEM_solver_1D_Poisson(left,right,h,basis_type_trial,basis_der_x_trial,basis_type_test,basis_der_x_test,basis_der_x_test_b,Gauss_type,s)
% basis_type_trial==101:1D linear
% N:��Ԫ����
% s:�����ĵ�����
%%%��������rightֵȡ�úܴ�Ļ������������������ֵ
%%%���һά�Ĳ��ɷ���
%%%�������Ա���������ﵽ�˿μ���Ҫ��
%%%�������ԣ���������ﵽ�˶��λ�����ʱ�����Ҫ��2021/4/27����������
%%
% clear
% clc
% close
% left=0;
% right=5;
% h=0.05;  %%%ʱ�䲽����ѡȡ����ʱ������д����Ϊ����뱣֤P������Ҷ˵㻹��1������ѡȡ�Ļ���MATLAB���ɵ��������ܱ�֤�Ҷ˵㻹��1
% ��ˣ�Ϊ�˱�֤�����ִ������⣬�ռ䲽����ʱȡ������1/12,1/24,1/48����  ������/2021/4/26��
% ����ǰ����Ĺ�����ʱ�ȷ�һ��;
% basis_type_trial=101;
% basis_type_test=101;
% basis_der_x_trial=1;
% basis_der_x_test=1;
% basis_der_x_test_b=0;
% Gauss_type=4;   %��˹�ڵ�ȡ�ĸ�
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
%%%�߽�ڵ�
boundary_type='DD';
boundarynodes=generate_boundarynodes(Pb,boundary_type);
%% ��װϵ������A
matrix_size=[max(size(Pb)) max(size(Pb))];
number_of_elements=size(T,2);
A=assemble_matrix_1D('function_c',Gauss_type,matrix_size,number_of_elements,P,T,Tb_trial,Tb_test,number_of_local_basis_fun_trial,number_of_local_basis_fun_test,basis_type_trial,basis_der_x_trial,basis_type_test,basis_der_x_test);

%%%%%��װ����b
%%
b=assemble_vector_1D('function_f',P,T,Tb,Nb,number_of_elements,number_of_local_basis_fun_test,Gauss_type,basis_type_test,basis_der_x_test_b);
%%%%�߽������Ĵ���
%%
[A0,b]=treat_Dirichlet_boundary('function_g',Pb,A,b,boundarynodes);
%%

% tic
% solution=pinv(full(A))*b;
solution=A0\b;

% solution=linsolve(full(A0),b);
% toc
% plot(P,solution,'*',P,P.*cos(P),'o')
%%%%�������
error=compute_Hs_error(P,T,Tb,s,solution,Gauss_type,'fx',number_of_local_basis_fun_test,basis_type_test);
end
% error=compute_max_FE_nodes(?);
% end