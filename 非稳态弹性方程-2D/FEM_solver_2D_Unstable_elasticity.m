 function [error,solution]=FEM_solver_2D_Unstable_elasticity(Start,End,dt,left,right,bottom,top,h1,h2,basis_type_trial,basis_type_test,basis_der_x_y_test_b,s)
 % basis_type_trial==201:2D linear
 % basis_type_trial==202:2D 二次元
N1=(right-left)/h1;
N2=(top-bottom)/h2;
 %%
[P,T]=generate_PT_2D(left,bottom,h1,h2,N1,N2);
[Pb,Tb]=generate_PbTb_2D(left,bottom,h1,h2,N1,N2,P,T,basis_type_trial);
%%
Gauss_type=8;                               %%%一维高斯节点的数目
Gauss_type_triangle=16;
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
boundaryedges=generate_boundaryedges(N1,N2,T,'D');
boundarynodes=generate_boundarynodes(N1,N2,basis_type_trial,'D');
%%
[xi,Ai]=generate_Gauss_reference_triangle(Gauss_type_triangle);     
[Ai_1D,Xi_1D]=generate_Gauss_reference_1D(Gauss_type);
A1=assemble_matrix_2D('function_lambda',matrix_size,number_of_elements,P,T,Tb_trial,Tb_test,number_of_local_basis_fun_trial,number_of_local_basis_fun_test,basis_type_trial,[1 0],basis_type_test,[1 0],xi,Ai);
A2=assemble_matrix_2D('function_nu',matrix_size,number_of_elements,P,T,Tb_trial,Tb_test,number_of_local_basis_fun_trial,number_of_local_basis_fun_test,basis_type_trial,[1 0],basis_type_test,[1 0],xi,Ai);
A3=assemble_matrix_2D('function_nu',matrix_size,number_of_elements,P,T,Tb_trial,Tb_test,number_of_local_basis_fun_trial,number_of_local_basis_fun_test,basis_type_trial,[0 1],basis_type_test,[0 1],xi,Ai);
A4=assemble_matrix_2D('function_lambda',matrix_size,number_of_elements,P,T,Tb_trial,Tb_test,number_of_local_basis_fun_trial,number_of_local_basis_fun_test,basis_type_trial,[0 1],basis_type_test,[1 0],xi,Ai);
A5=assemble_matrix_2D('function_nu',matrix_size,number_of_elements,P,T,Tb_trial,Tb_test,number_of_local_basis_fun_trial,number_of_local_basis_fun_test,basis_type_trial,[1 0],basis_type_test,[0 1],xi,Ai);
A6=assemble_matrix_2D('function_lambda',matrix_size,number_of_elements,P,T,Tb_trial,Tb_test,number_of_local_basis_fun_trial,number_of_local_basis_fun_test,basis_type_trial,[1 0],basis_type_test,[0 1],xi,Ai);
A7=assemble_matrix_2D('function_nu',matrix_size,number_of_elements,P,T,Tb_trial,Tb_test,number_of_local_basis_fun_trial,number_of_local_basis_fun_test,basis_type_trial,[0 1],basis_type_test,[1 0],xi,Ai);
A8=assemble_matrix_2D('function_lambda',matrix_size,number_of_elements,P,T,Tb_trial,Tb_test,number_of_local_basis_fun_trial,number_of_local_basis_fun_test,basis_type_trial,[0 1],basis_type_test,[0 1],xi,Ai);

A=[A1+2*A2+A3 A4+A5; A6+A7 A8+2*A3+A2];

Me=assemble_matrix_2D('function_one',matrix_size,number_of_elements,P,T,Tb_trial,Tb_test,number_of_local_basis_fun_trial,number_of_local_basis_fun_test,basis_type_trial,[0 0],basis_type_test,[0 0],xi,Ai);
O4=sparse(Nb,Nb);

M=[Me,O4;O4,Me];

Atilde=M/dt/dt+A/4;
Afixed_1=2*M/dt/dt-A/2;
Afixed_2=  M/dt/dt+A/4;
clear M A;
%%
dt=1;
[X_1,X_0]=generate_initial_vector(Pb,dt);
Xm=X_1;  
Xm_1=X_0;
% an1=fx1(Pb(1,:),Pb(2,:),dt,[0 0])';  %%analytical solution
% an2=fx2(Pb(1,:),Pb(2,:),dt,[0 0])';
% an=[an1;an2];
% plot(X_1,an,'o')
%%
number_of_time_step=(End-Start)/dt;

for m=1:number_of_time_step-1
    tm=m*dt;
    
    tmp1=(m+1)*dt;
    
    b1m=assemble_vector_2D_t('function_f1',tm,P,T,Tb,matrix_size(2),number_of_elements,number_of_local_basis_fun_test,basis_type_test,basis_der_x_y_test_b,xi,Ai);
    
    b2m=assemble_vector_2D_t('function_f2',tm,P,T,Tb,matrix_size(2),number_of_elements,number_of_local_basis_fun_test,basis_type_test,basis_der_x_y_test_b,xi,Ai);
    
    bm=[b1m;b2m];

    btilde=bm+Afixed_1*Xm-Afixed_2*Xm_1;
    
    [Atilde,btilde]=treat_Dirichlet_boundary_t('function_b1','function_b2',tmp1,Pb,Atilde,btilde,boundarynodes);
    
    Xmp1=Atilde\btilde;
    
    Xm_1=Xm;
    
    Xm=Xmp1;
end
solution=Xm;
% v1=treat_Neumann_boundary_2D('function_p1',P,T,Pb,Tb,Ai_1D,Xi_1D,boundaryedges,number_of_local_basis_fun_test,basis_type_test);
% v2=treat_Neumann_boundary_2D('function_p2',P,T,Pb,Tb,Ai_1D,Xi_1D,boundaryedges,number_of_local_basis_fun_test,basis_type_test);
% b=b+[v1;v2];
% [w,R]=treat_Robin_boundary_2D('function_cq','function_cr',P,T,Pb,Tb,Ai_1D,Xi_1D,boundaryedges,number_of_local_basis_fun_trial,number_of_local_basis_fun_test,basis_type_trial,basis_type_test);
% A=A+R;
% b=b+w;
%%%Dirichlet边界条件一定要最后处理
% [A,b]=treat_Dirichlet_boundary('function_b1','function_b2',Pb,A,b,boundarynodes);
% Neumann 边界条件

%%
an1=fx1(Pb(1,:),Pb(2,:),End,[0 0])';  %%analytical solution
an2=fx2(Pb(1,:),Pb(2,:),End,[0 0])';
an=[an1;an2];
%%
plot(solution,an,'*')
% error=max(abs(solution-an));
%% Compute the error
%  给出高斯节点上的无穷范数
error_inf=max(abs(an-solution));
error=compute_Hs_error_2D(P,T,Tb,s,solution,'fx1','fx2',End,number_of_local_basis_fun_test,basis_type_test,xi,Ai);
end