function result=Gauss_int_line_trial_test(function_cr,Gauss_weights,Gauss_nodes,vertices,basis_type_trial,basis_type_test,basis_index_trial,basis_index_test,basis_der_x,basis_der_y)


Gpn=length(Gauss_weights);

int_value=0;

for k=1:Gpn
   int_value=int_value+Gauss_weights(k)*feval(function_cr,Gauss_nodes(1,k),Gauss_nodes(2,k))...
       *local_basis_2D(Gauss_nodes(1,k),Gauss_nodes(2,k),vertices,basis_type_trial,basis_index_trial,basis_der_x,basis_der_y)...
   *local_basis_2D(Gauss_nodes(1,k),Gauss_nodes(2,k),vertices,basis_type_test,basis_index_test,basis_der_x,basis_der_y);
end

result=int_value;
end