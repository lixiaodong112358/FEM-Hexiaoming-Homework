function result=Gauss_int_error_2D(Gauss_weight,Gauss_nodes,End,s,uh_local_vec1,analytical_solution1,NLB,vertices,basis_type)
%%% End 最后的时刻
Gpn=size(Gauss_nodes,2);%%高斯积分点
int_value=0;
for k=1:Gpn
   int_value=int_value+...
       Gauss_weight(k)*(feval(analytical_solution1,Gauss_nodes(k,1),Gauss_nodes(k,2),End,[s(1) 0])-local_FE_function_2D(Gauss_nodes(k,1),Gauss_nodes(k,2),uh_local_vec1,NLB,vertices,basis_type,[s(1) 0]))^2+...
       Gauss_weight(k)*(feval(analytical_solution1,Gauss_nodes(k,1),Gauss_nodes(k,2),End,[0 s(2)])-local_FE_function_2D(Gauss_nodes(k,1),Gauss_nodes(k,2),uh_local_vec1,NLB,vertices,basis_type,[0 s(2)]))^2;
end
result=int_value;
end