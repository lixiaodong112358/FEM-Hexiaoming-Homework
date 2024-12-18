 function [error_u,error_p,solution]=FEM_solver_2D_Unsteady_stokes(theta,Start,End,dt,left,right,bottom,top,h1,h2,basis_type_u,basis_der_x_y_trial,basis_type_p,basis_der_x_y_test,basis_der_x_y_test_b,s)
 % basis_type_trial==201:2D linear
 % basis_type_trial==202:2D 二次元
N1=(right-left)/h1;
N2=(top-bottom)/h2;

 %%
 basis_type_trial=basis_type_u;
[P,T]=generate_PT_2D(left,bottom,h1,h2,N1,N2);
[Pb,Tb]=generate_PbTb_2D(left,bottom,h1,h2,N1,N2,P,T,basis_type_trial);
%%
Gauss_type=8;                               %%%一维高斯节点的数目
Gauss_type_triangle=9;
matrix_size=[size(Pb,2),size(Pb,2)];
number_of_elements=size(T,2);
N=size(P,2);
Nb=size(Pb,2);
Nbp=N;
%%
if basis_type_u==201
  Pb_u=P;
  Tb_u=T;
  number_of_local_basis_fun_u=3;
elseif basis_type_u==202
  Pb_u=Pb;
  Tb_u=Tb;
  number_of_local_basis_fun_u=6;
end
if basis_type_p==201
  Pb_p=P;
  Tb_p=T;
  number_of_local_basis_fun_p=3;
elseif basis_type_p==202
  Pb_p=Pb;
  Tb_p=Tb;
  number_of_local_basis_fun_p=6;
end
%%
% boundaryedges=generate_boundaryedges(N1,N2,Tb,'D');
% boundarynodes=generate_boundarynodes(N1,N2,basis_type_trial,'D');
boundaryedges=generate_boundaryedges(N1,N2,T,'D');
boundarynodes=generate_boundarynodes(N1,N2,basis_type_trial,'D');
%%
[xi,Ai]=generate_Gauss_reference_triangle(Gauss_type_triangle);     
[Ai_1D,Xi_1D]=generate_Gauss_reference_1D(Gauss_type);
A1=assemble_matrix_2D('function_v',matrix_size,number_of_elements,P,T,Tb_u,Tb_u,number_of_local_basis_fun_u,number_of_local_basis_fun_u,basis_type_u,[1 0],basis_type_u,[1 0],xi,Ai);
A2=assemble_matrix_2D('function_v',matrix_size,number_of_elements,P,T,Tb_u,Tb_u,number_of_local_basis_fun_u,number_of_local_basis_fun_u,basis_type_u,[0 1],basis_type_u,[0 1],xi,Ai);
A3=assemble_matrix_2D('function_v',matrix_size,number_of_elements,P,T,Tb_u,Tb_u,number_of_local_basis_fun_u,number_of_local_basis_fun_u,basis_type_u,[1 0],basis_type_u,[0 1],xi,Ai);
A5=assemble_matrix_2D('function_negative_one',[Nb Nbp],number_of_elements,P,T,Tb_p,Tb_u,number_of_local_basis_fun_p,number_of_local_basis_fun_u,basis_type_p,[0 0],basis_type_u,[1 0],xi,Ai);
A6=assemble_matrix_2D('function_negative_one',[Nb Nbp],number_of_elements,P,T,Tb_p,Tb_u,number_of_local_basis_fun_p,number_of_local_basis_fun_u,basis_type_p,[0 0],basis_type_u,[0 1],xi,Ai);
O1=sparse(Nbp,Nbp);
A=[2*A1+A2 A3 A5; A3' 2*A2+A1 A6;A5' A6' O1];

Me=assemble_matrix_2D('function_Me',matrix_size,number_of_elements,P,T,Tb_u,Tb_u,number_of_local_basis_fun_u,number_of_local_basis_fun_u,basis_type_u,[0 0],basis_type_u,[0 0],xi,Ai);
O2=sparse(Nb,Nbp);
O3=sparse(Nb,Nb);
M=[Me O3 O2;O3 Me O2;O2' O2' O1];

Atilde=M/dt+theta*A;
Afixed=M/dt-(1-theta)*A;
clear M A;

X_0=generate_initial_vector('function_initial_u1','function_initial_u2','function_initial_p',Pb,P);
X_old=X_0;  
number_of_time_step=(End-Start)/dt;
b0=zeros(Nbp,1);
for m=0:number_of_time_step-1
    tm=m*dt;
    
    tmp1=(m+1)*dt;
    
    b1m=assemble_vector_2D_t('function_f1',tm  ,P,T,Tb,matrix_size(2),number_of_elements,number_of_local_basis_fun_u,basis_type_u,basis_der_x_y_test_b,xi,Ai);
    b2m=assemble_vector_2D_t('function_f2',tm  ,P,T,Tb,matrix_size(2),number_of_elements,number_of_local_basis_fun_u,basis_type_u,basis_der_x_y_test_b,xi,Ai);
    bm=[b1m;b2m;b0];
    b1mp1=assemble_vector_2D_t('function_f1',tmp1,P,T,Tb,matrix_size(2),number_of_elements,number_of_local_basis_fun_u,basis_type_u,basis_der_x_y_test_b,xi,Ai);
    b2mp1=assemble_vector_2D_t('function_f2',tmp1,P,T,Tb,matrix_size(2),number_of_elements,number_of_local_basis_fun_u,basis_type_u,basis_der_x_y_test_b,xi,Ai);
    bmp1=[b1mp1;b2mp1;b0];
    
    btilde=theta*bmp1+(1-theta)*bm+Afixed*X_old;
    
    [Atilde,btilde]=fix_p('fxp',tmp1,Atilde,btilde,P,Nb);
    
    [Atilde,btilde]=treat_Dirichlet_boundary_t('function_b1','function_b2',tmp1,Pb,Atilde,btilde,boundarynodes);
     
    solution=Atilde\btilde;
    

    
%     solution-an
%     result=[solution,result];
    
    X_old=solution;
end
% an=fx(Pb(1,:),Pb(2,:),tmp1,[0 0])';
% b1=assemble_vector_2D('function_f1',P,T,Tb_u,matrix_size(2),number_of_elements,number_of_local_basis_fun_u,basis_type_u,basis_der_x_y_test_b,xi,Ai);
% b2=assemble_vector_2D('function_f2',P,T,Tb_u,matrix_size(2),number_of_elements,number_of_local_basis_fun_u,basis_type_u,basis_der_x_y_test_b,xi,Ai);
% b0=zeros(Nbp,1);
% b=[b1;b2;b0];
% v1=treat_Neumann_boundary_2D('function_p1',P,T,Pb,Tb,Ai_1D,Xi_1D,boundaryedges,number_of_local_basis_fun_test,basis_type_test);
% v2=treat_Neumann_boundary_2D('function_p2',P,T,Pb,Tb,Ai_1D,Xi_1D,boundaryedges,number_of_local_basis_fun_test,basis_type_test);
% b=b+[v1;v2;b0];
% [w,R]=treat_Robin_boundary_2D('function_cq','function_cr',P,T,Pb,Tb,Ai_1D,Xi_1D,boundaryedges,number_of_local_basis_fun_trial,number_of_local_basis_fun_test,basis_type_trial,basis_type_test);
% A=A+R;
% b=b+w;
%%%Dirichlet边界条件一定要最后处理
% [A,b]=fix_p('fxp',A,b,P,Nb);
% [A,b]=treat_Dirichlet_boundary('function_b1','function_b2',Pb,A,b,boundarynodes);
% Neumann 边界条件

%%
% solution=A\b;
%%最后时刻的解析解
an1=fx1(Pb(1,:),Pb(2,:),End,[0 0])';  %%analytical solution
an2=fx2(Pb(1,:),Pb(2,:),End,[0 0])';
an3=fxp( P(1,:), P(2,:),End,[0 0])';
an=[an1;an2;an3];
%%
figure
plot(solution(1:2*Nb),an(1:2*Nb),'*')
title('位移值解析解和数值解对比')
% error=max(abs(solution-an));
%%
figure
plot(solution(2*Nb+1:end),an3,'o')
title('压力值P解析解和数值解对比')
%% Compute the error 计算最后时刻的误差
%  给出高斯节点上的无穷范数
% error_inf=max(abs(an-solution));
error_u=compute_Hs_error_2D(P,T,Tb_u,End,s,solution(1:2*Nb),'fx1','fx2',number_of_local_basis_fun_u,basis_type_u,xi,Ai);
error_p=compute_Hs_error_2D_p(P,T,Tb_p,End,s,solution(2*Nb+1:end),'fxp',number_of_local_basis_fun_p,basis_type_p,xi,Ai);
end