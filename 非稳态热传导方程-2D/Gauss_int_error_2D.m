function result=Gauss_int_error_2D(t,Gauss_weight,Gauss_nodes,s,uh_local_vec,analytical_solution,NLB,vertices,basis_type)

Gpn=size(Gauss_nodes,2);%%高斯积分点
int_value=0;
for k=1:Gpn
   int_value=int_value+Gauss_weight(k)*(feval(analytical_solution,Gauss_nodes(k,1),Gauss_nodes(k,2),t,s)-local_FE_function_2D(Gauss_nodes(k,1),Gauss_nodes(k,2),uh_local_vec,NLB,vertices,basis_type,s))^2;
end
result=int_value;
end