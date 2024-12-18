function [result]=assemble_vector_1D(coe_fun,P,T,Tb,Nb,number_of_elements,NLB,Gauss_type,basis_type_test,basis_der_x_test)
b=zeros(Nb,1);
for n=1:number_of_elements
    vertices=P(:,T(:,n));
    [Gauss_weight,Gauss_nodes]=generate_Gauss_local_1D(vertices,Gauss_type);
   for beta=1:NLB
       basis_index_test=beta;
       r=Gauss_quad_1D_test(coe_fun,Gauss_weight,Gauss_nodes,vertices,basis_type_test,basis_index_test,basis_der_x_test);
       b(Tb(beta,n),1)=b(Tb(beta,n),1)+r;
   end
end
result=b;
end