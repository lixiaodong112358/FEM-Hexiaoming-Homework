function A=assemble_matrix_1D(coe_fun,Gauss_type,matrix_size,number_of_elements,P,T,Tb_trial,Tb_test,number_of_local_basis_fun_trial,number_of_local_basis_fun_test,basis_type_trial,basis_der_x_trial,basis_type_test,basis_der_x_test)
%%%%%组装总刚度矩阵
%%%%%coe_fun:系数函数的函数名
%%%%%2021/4/19 此函数初步测试完毕 值得庆祝--李晓东
A=sparse(matrix_size(1),matrix_size(2));

for n=1:number_of_elements
    
    vertices=P(:,T(:,n));
     
    [Gauss_weight,Gauss_nodes]=generate_Gauss_local_1D(vertices,Gauss_type);
    

    for alpha=1:number_of_local_basis_fun_trial
        for beta=1:number_of_local_basis_fun_test
            basis_index_trial=alpha;
            basis_index_test=beta;
            int_value=Gauss_quad_1D_trial_test(coe_fun,Gauss_weight,Gauss_nodes,vertices,basis_type_trial,basis_index_trial,basis_der_x_trial,basis_type_test,basis_index_test,basis_der_x_test);
%             S(beta,alpha)=int_value
            i=Tb_test(beta,n);
            j=Tb_trial(alpha,n);
            A(i,j)=A(i,j)+int_value;
        end
    end
        

end