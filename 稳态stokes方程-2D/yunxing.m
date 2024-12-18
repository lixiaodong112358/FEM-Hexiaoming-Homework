%%%%本程序用于测试每一个函数，2021/5/14
clc
clear
close

left=0;
right=1;
bottom=0;
top=1;
N1=2;
N2=2;
basis_type=201;
basis_type_trial=201;
basis_type_test=201;
%%
h1=(right-left)/N1;
h2=(top-bottom)/N2;
% mesh=generate_mesh(N1,N2,basis_type_trial);
[P,T]=generate_PT_2D(left,bottom,h1,h2,N1,N2);
[Pb,Tb]=generate_PbTb_2D(left,bottom,h1,h2,N1,N2,P,T,basis_type_trial);

coe_fun='function_c';
Gauss_type=4;
matrix_size=[size(Pb,2),size(Pb,2)];
number_of_elements=size(T,2);

if basis_type_trial==201
  Pb_trial=P;
  Tb_trial=T;
  number_of_local_basis_fun_trial=3;
elseif basis_type_trial==202
  Pb_trial=Pb;
  Tb_trial=Tb;
  number_of_local_basis_fun_trial=6;
end
if basis_type_test==201
  Pb_test=P;
  Tb_test=T;
  number_of_local_basis_fun_test=3;
elseif basis_type_test==202
  Pb_test=Pb;
  Tb_test=Tb;
  number_of_local_basis_fun_test=6;
end
basis_der_x_y_trial=[1 0;0 1];
basis_der_x_y_test=[1 0;0 1];
basis_der_x_y_test_b=[0 0];
%%
boundaryedges=generate_boundaryedges(N1,N2,Tb,'D');
boundarynodes=generate_boundarynodes(N1,N2,basis_type_trial,'D');
%%
[Gauss_point_reference_triangle,Gauss_coefficient_reference_triangle]=generate_Gauss_reference_triangle(Gauss_type);
A=assemble_matrix_2D('function_c',matrix_size,number_of_elements,P,T,Tb_trial,Tb_test,number_of_local_basis_fun_trial,number_of_local_basis_fun_test,basis_type_trial,basis_der_x_y_trial,basis_type_test,basis_der_x_y_test,Gauss_point_reference_triangle,Gauss_coefficient_reference_triangle);
% full(A)
b=assemble_vector_2D('function_f',P,T,Tb,matrix_size(2),number_of_elements,number_of_local_basis_fun_test,basis_type_test,basis_der_x_y_test_b,Gauss_point_reference_triangle,Gauss_coefficient_reference_triangle);
[A,b]=treat_Dirichlet_boundary('function_b',Pb,A,b,boundarynodes);
%%
solution=A\b