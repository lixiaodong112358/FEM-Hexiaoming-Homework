function A=assemble_matrix_2D_coe(coe_fun,uh_vec,basis_type_coe,basis_der_x_coe,basis_der_y_coe,matrix_size,number_of_elements,P,T,Tb_trial,Tb_test,number_of_local_basis_fun_trial,number_of_local_basis_fun_test,basis_type_trial,basis_der_x_y_trial,basis_type_test,basis_der_x_y_test,Gauss_point_reference_triangle,Gauss_coefficient_reference_triangle)
%%%%% 组装总刚度矩阵
%%%%% NS方程非线性项需要的组装器
%%%%% coe_fun:系数函数的函数名
%%%%% 2021/7/19 
%%
A=sparse(matrix_size(1),matrix_size(2));
NLB=number_of_local_basis_fun_trial;
for n=1:number_of_elements

    vertices=P(:,T(:,n));
      
    [Gauss_weight,Gauss_nodes]=generate_Gauss_local_triangle(Gauss_coefficient_reference_triangle,Gauss_point_reference_triangle,vertices);

    uh_local_vec=uh_vec(Tb_trial(:,n));
    
    for alpha=1:number_of_local_basis_fun_trial
        for beta=1:number_of_local_basis_fun_test
            basis_index_trial=alpha;
            basis_index_test=beta;
            int_value=Gauss_quad_2D_trial_test_coe(coe_fun,uh_local_vec,basis_type_coe,NLB,basis_der_x_coe,basis_der_y_coe,Gauss_weight,Gauss_nodes,vertices,basis_type_trial,basis_index_trial,basis_type_test,basis_index_test,basis_der_x_y_test,basis_der_x_y_trial);
%             S(beta,alpha)=int_value
            i=Tb_test(beta,n);
            j=Tb_trial(alpha,n);
            A(i,j)=A(i,j)+int_value;
        end
    end
        

end