 function [error,solution]=FEM_solver_2D_Poisson(left,right,bottom,top,h1,h2,basis_type_trial,basis_der_x_y_trial,basis_type_test,basis_der_x_y_test,basis_der_x_y_test_b,s)
 % basis_type_trial==201:2D linear
 % basis_type_trial==202:2D ถþดฮิช
N1=(right-left)/h1;
N2=(top-bottom)/h2;
 %%
[P,T]=generate_PT_2D(left,bottom,h1,h2,N1,N2);
[Pb,Tb]=generate_PbTb_2D(left,bottom,h1,h2,N1,N2,P,T,basis_type_trial);
%%
Gauss_type=4;
matrix_size=[size(Pb,2),size(Pb,2)];
number_of_elements=size(T,2);
N=size(P,2);
Nb=size(Pb,2);
%%
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
%%
% boundaryedges=generate_boundaryedges(N1,N2,Tb,'D');
% boundarynodes=generate_boundarynodes(N1,N2,basis_type_trial,'D');
boundaryedges=generate_boundaryedges(N1,N2,Tb,'D');
boundarynodes=generate_boundarynodes(N1,N2,basis_type_trial,'D');
%%
[xi,Ai]=generate_Gauss_reference_triangle(Gauss_type);
A1=assemble_matrix_2D('function_c',matrix_size,number_of_elements,P,T,Tb_trial,Tb_test,number_of_local_basis_fun_trial,number_of_local_basis_fun_test,basis_type_trial,[1 0],basis_type_test,[1 0],xi,Ai);
A2=assemble_matrix_2D('function_c',matrix_size,number_of_elements,P,T,Tb_trial,Tb_test,number_of_local_basis_fun_trial,number_of_local_basis_fun_test,basis_type_trial,[0 1],basis_type_test,[0 1],xi,Ai);
A=A1+A2;
b=assemble_vector_2D('function_f',P,T,Tb,matrix_size(2),number_of_elements,number_of_local_basis_fun_test,basis_type_test,basis_der_x_y_test_b,xi,Ai);
v=treat_Neumann_boundary_2D('function_cp',P,T,Pb,Tb,boundaryedges,number_of_local_basis_fun_test,basis_type_test);
b=b+v;
[A,b]=treat_Dirichlet_boundary('function_b',Pb,A,b,boundarynodes);
%%
solution=A\b;
an=fx(Pb(1,:),Pb(2,:),[0 0])';  %%analytical solution
error=max(abs(solution-an))
plot(an,solution,'*')
%% Compute the error
error_1=compute_Hs_error_2D(P,T,Tb,[1 1],solution,'fx',number_of_local_basis_fun_test,basis_type_test,xi,Ai);
error_0=compute_Hs_error_2D(P,T,Tb,[0 0],solution,'fx',number_of_local_basis_fun_test,basis_type_test,xi,Ai);
end