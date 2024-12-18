function int_value=Gauss_quad_2D_trial_test_coe(coe_fun,uh_local_vec,basis_type_coe,NLB,basis_der_x_coe,basis_der_y_coe,Gauss_weight,Gauss_nodes,vertices,basis_type_trial,basis_index_trial,basis_type_test,basis_index_test,basis_der_x_y_test,basis_der_x_y_trial)
%%%处理非线性项需要的高斯积分
r=basis_der_x_y_trial(1);
s=basis_der_x_y_trial(2);
p=basis_der_x_y_test(1);
q=basis_der_x_y_test(2);
Gpn=length(Gauss_weight);%%高斯积分点
int_value=0;
for k=1:Gpn
   int_value=int_value+Gauss_weight(k)*feval(coe_fun,Gauss_nodes(k,1),Gauss_nodes(k,2),uh_local_vec,NLB,vertices,basis_type_coe,[basis_der_x_coe,basis_der_y_coe])*...
       local_basis_2D(Gauss_nodes(k,1),Gauss_nodes(k,2),vertices,basis_type_trial,basis_index_trial,r,s)*...
       local_basis_2D(Gauss_nodes(k,1),Gauss_nodes(k,2),vertices,basis_type_test,basis_index_test,p,q);
end
end