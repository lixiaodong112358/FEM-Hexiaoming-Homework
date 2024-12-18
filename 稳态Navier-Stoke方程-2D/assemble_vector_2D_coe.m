function [result]=assemble_vector_2D_coe(coe_fun,u1h_vec,u2h_vec,basis_type_coe1,basis_der_x_coe1,basis_der_y_coe1,basis_type_coe2,basis_der_x_coe2,basis_der_y_coe2,P,T,Tb,Nb,number_of_elements,number_of_local_basis_fun_test,basis_type_test,basis_der_x_y_test_b,Gauss_point_reference_triangle,Gauss_coefficient_reference_triangle)
b=zeros(Nb,1);
NLB=number_of_local_basis_fun_test;
for n=1:number_of_elements
    vertices=P(:,T(:,n));
    
    [Gauss_weight,Gauss_nodes]=generate_Gauss_local_triangle(Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,vertices);
    
    u1h_local_vec=u1h_vec(Tb(:,n));
    
    u2h_local_vec=u2h_vec(Tb(:,n));
    
   for beta=1:number_of_local_basis_fun_test
       
       basis_index_test=beta;
       
       r=Gauss_quad_2D_b_coe(coe_fun,NLB,u1h_local_vec,u2h_local_vec,basis_type_coe1,basis_der_x_coe1,basis_der_y_coe1,basis_type_coe2,basis_der_x_coe2,basis_der_y_coe2,Gauss_weight,Gauss_nodes,vertices,basis_type_test,basis_index_test,basis_der_x_y_test_b);
       
       b(Tb(beta,n),1)=b(Tb(beta,n),1)+r;
   end
end
result=b;
end