function [result]=Gauss_quad_2D_b_coe(coe_fun,NLB,u1h_local_vec,u2h_local_vec,basis_type_coe1,basis_der_x_coe1,basis_der_y_coe1,basis_type_coe2,basis_der_x_coe2,basis_der_y_coe2,Gauss_weight,Gauss_nodes,vertices,basis_type_test,basis_index_test,basis_der_x_y_test_b)


Gpn=length(Gauss_nodes);%%高斯积分点
int_value=0;
for k=1:Gpn
   int_value=int_value+Gauss_weight(k)*feval(coe_fun,Gauss_nodes(k,1),Gauss_nodes(k,2),u1h_local_vec,NLB,vertices,basis_type_coe1,[basis_der_x_coe1,basis_der_y_coe1])*...
       feval(coe_fun,Gauss_nodes(k,1),Gauss_nodes(k,2),u2h_local_vec,NLB,vertices,basis_type_coe2,[basis_der_x_coe2,basis_der_y_coe2])*...
   local_basis_2D(Gauss_nodes(k,1),Gauss_nodes(k,2),vertices,basis_type_test,basis_index_test,basis_der_x_y_test_b(1),basis_der_x_y_test_b(2));
end
result=int_value;
end