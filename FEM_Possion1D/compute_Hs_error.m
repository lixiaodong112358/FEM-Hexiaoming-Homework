function error=compute_Hs_error(P,T,Tb,s,solution_vector,Gauss_type,analytical_solution,NLB,basis_type)
%%%s�ǵ�����
%%%solution_vec����ֵ������
%%%analytical_solution��������ĺ�����
number_of_elements=size(T,2);
error=0;
for n=1:number_of_elements
    
    vertices=P(:,T(:,n));
     
    [Gauss_weight,Gauss_nodes]=generate_Gauss_local_1D(vertices,Gauss_type);
    
    uh_local_vec=solution_vector(Tb(:,n));
    
    int_value=Gauss_int_error_1D(Gauss_weight,Gauss_nodes,s,uh_local_vec,analytical_solution,NLB,vertices,basis_type);
    
    error=error+int_value;
     
end
error=sqrt(error);
end