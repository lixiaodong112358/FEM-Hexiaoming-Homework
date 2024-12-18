 function [error,solution]=FEM_solver_2D_heat(theta,Start,End,dt,left,right,bottom,top,h1,h2,basis_type_trial,basis_der_x_y_trial,basis_type_test,basis_der_x_y_test,basis_der_x_y_test_b,s)
 % basis_type_trial==201:2D linear
 % basis_type_trial==202:2D 二次元
 % 对于一个矩形区域求解非稳态热传导问题  
 % Start 开始时刻
 % End  :结束时刻
N1=(right-left)/h1;
N2=(top-bottom)/h2;
 %%
[P,T]=generate_PT_2D(left,bottom,h1,h2,N1,N2);
[Pb,Tb]=generate_PbTb_2D(left,bottom,h1,h2,N1,N2,P,T,basis_type_trial);
%%
Gauss_type=8;                               %%%一维高斯节点的数目
Gauss_type_triangle=9;
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
[xi,Ai]=generate_Gauss_reference_triangle(Gauss_type_triangle);     
[Ai_1D,Xi_1D]=generate_Gauss_reference_1D(Gauss_type);% %%%处理Neumann边界条件和Robin 边界条件所用
M=assemble_matrix_2D('function_M',matrix_size,number_of_elements,P,T,Tb_trial,Tb_test,number_of_local_basis_fun_trial,number_of_local_basis_fun_test,basis_type_trial,[0 0],basis_type_test,[0 0],xi,Ai);
A1=assemble_matrix_2D('function_c',matrix_size,number_of_elements,P,T,Tb_trial,Tb_test,number_of_local_basis_fun_trial,number_of_local_basis_fun_test,basis_type_trial,[1 0],basis_type_test,[1 0],xi,Ai);
A2=assemble_matrix_2D('function_c',matrix_size,number_of_elements,P,T,Tb_trial,Tb_test,number_of_local_basis_fun_trial,number_of_local_basis_fun_test,basis_type_trial,[0 1],basis_type_test,[0 1],xi,Ai);
A=A1+A2;
% theta=0.8;                  %%Grank-Nicolson scheme
Atilde=M/dt+theta*A;
Afixed=M/dt-(1-theta)*A;
clear M A;
X_0=generate_initial_vector('function_initial',Pb);
X_old=X_0;  
number_of_time_step=(End-Start)/dt;
%%%%   时间迭代：
% result=[];
for m=0:number_of_time_step-1
    tm=m*dt;
    
    tmp1=(m+1)*dt;
    
    bm=assemble_vector_2D_t('function_f',tm,P,T,Tb,matrix_size(2),number_of_elements,number_of_local_basis_fun_test,basis_type_test,basis_der_x_y_test_b,xi,Ai);
    
    bmp1=assemble_vector_2D_t('function_f',tmp1,P,T,Tb,matrix_size(2),number_of_elements,number_of_local_basis_fun_test,basis_type_test,basis_der_x_y_test_b,xi,Ai);
    
    btilde=theta*bmp1+(1-theta)*bm+Afixed*X_old;
    
    [Atilde,btilde]=treat_Dirichlet_boundary_t('function_b',tmp1,Pb,Atilde,btilde,boundarynodes);
    
    v=treat_Neumann_boundary_2D('function_cp',tmp1,P,T,Pb,Tb,Ai_1D,Xi_1D,boundaryedges,number_of_local_basis_fun_test,basis_type_test);
    
    btilde=btilde+v;
    
    [w,R]=treat_Robin_boundary_2D('function_cq','function_cr',P,T,Pb,Tb,Ai_1D,Xi_1D,boundaryedges,number_of_local_basis_fun_trial,number_of_local_basis_fun_test,basis_type_trial,basis_type_test);
    
    Atilde=Atilde+R;
    
    btilde=btilde+w;
    
    solution=Atilde\btilde;
    
    an=fx(Pb(1,:),Pb(2,:),tmp1,[0 0])';
    
%     solution-an
%     result=[solution,result];
    
    X_old=solution;
end
% b=assemble_vector_2D('function_f',P,T,Tb,matrix_size(2),number_of_elements,number_of_local_basis_fun_test,basis_type_test,basis_der_x_y_test_b,xi,Ai);
% v=treat_Neumann_boundary_2D('function_cp',P,T,Pb,Tb,Ai_1D,Xi_1D,boundaryedges,number_of_local_basis_fun_test,basis_type_test);
% b=b+v;
% [w,R]=treat_Robin_boundary_2D('function_cq','function_cr',P,T,Pb,Tb,Ai_1D,Xi_1D,boundaryedges,number_of_local_basis_fun_trial,number_of_local_basis_fun_test,basis_type_trial,basis_type_test);
% A=A+R;
% b=b+w;
% [A,b]=treat_Dirichlet_boundary_t('function_b',Pb,A,b,boundarynodes);
% Neumann 边界条件

%%
% solution=A\b;
an0=fx(Pb(1,:),Pb(2,:),0,[0 0])'; 
an=fx(Pb(1,:),Pb(2,:),tmp1,[0 0])';  %%analytical solution
%%
plot(solution,an,'*')
figure
plot(1:Nb,solution,'o:',1:Nb,an,'r--')
% figure
% plot(1:Nb,X_0,'o:',1:Nb,an0,'p--')
% error=max(abs(solution-an));
%% Compute the error
error=compute_Hs_error_2D(tmp1,P,T,Tb,s,solution,'fx',number_of_local_basis_fun_test,basis_type_test,xi,Ai);
end