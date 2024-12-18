function [result]=assemble_vector_2D_t(coe_fun,t,P,T,Tb,Nb,number_of_elements,number_of_local_basis_fun_test,basis_type_test,basis_der_x_y_test_b,Gauss_point_reference_triangle,Gauss_coefficient_reference_triangle)
b=zeros(Nb,1);

for n=1:number_of_elements
    vertices=P(:,T(:,n));
    [Gauss_weight,Gauss_nodes]=generate_Gauss_local_triangle(Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,vertices);
   for beta=1:number_of_local_basis_fun_test
       basis_index_test=beta;
       r=Gauss_quad_2D_test(coe_fun,t,Gauss_weight,Gauss_nodes,vertices,basis_type_test,basis_index_test,basis_der_x_y_test_b);
       b(Tb(beta,n),1)=b(Tb(beta,n),1)+r;
   end
end
result=b;
end