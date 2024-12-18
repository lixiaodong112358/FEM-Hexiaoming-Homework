function error=compute_Hs_error_2D(P,T,Tb_trial,s,solution_vector,analytical_solution,NLB,basis_type,Gauss_point_reference_triangle,Gauss_coefficient_reference_triangle)
%%%s�ǵ�����
%%%solution_vec����ֵ������
%%%analytical_solution��������ĺ�����
number_of_elements=size(T,2);
error=0;
for n=1:number_of_elements
    
    vertices=P(:,T(:,n));
     
    [Gauss_weight,Gauss_nodes]=generate_Gauss_local_triangle(Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,vertices);
    
    uh_local_vec=solution_vector(Tb_trial(:,n));
    
    if s(1)==0&&s(2)==0
    
    int_value=Gauss_int_error_2D(Gauss_weight,Gauss_nodes,s,uh_local_vec,analytical_solution,NLB,vertices,basis_type);
    
    elseif s(1)==1&&s(2)==1
        
    int_value=Gauss_int_error_2D(Gauss_weight,Gauss_nodes,[0 s(2)],uh_local_vec,analytical_solution,NLB,vertices,basis_type)+Gauss_int_error_2D(Gauss_weight,Gauss_nodes,[s(1) 0],uh_local_vec,analytical_solution,NLB,vertices,basis_type);
    
    else
        warning='�����������������������'
    end
    
    error=error+int_value;
     
end
error=sqrt(error);
end