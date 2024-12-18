function int_value=Gauss_int_line_test(Neumann_fun,t,Gauss_weights,Gauss_nodes,vertices,basis_type,basis_index,basis_der_x,basis_der_y)

Gpn=length(Gauss_weights);

int_value=0;

for k=1:Gpn
   int_value=int_value+Gauss_weights(k)*feval(Neumann_fun,Gauss_nodes(1,k),Gauss_nodes(2,k),t)...
       *local_basis_2D(Gauss_nodes(1,k),Gauss_nodes(2,k),vertices,basis_type,basis_index,basis_der_x,basis_der_y);
end



end