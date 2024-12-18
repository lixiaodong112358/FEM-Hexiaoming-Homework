function A=assemble_matrix_2D_CVFEM(coe_fun,matrix_size,number_of_elements,P,T,Tb_trial,Tb_test,number_of_local_basis_fun_trial,number_of_local_basis_fun_test,basis_type_trial,basis_der_x_y_trial,basis_type_test,basis_der_x_y_test,Gauss_point_reference_triangle,Gauss_coefficient_reference_triangle)


A=sparse(matrix_size(1),matrix_size(2));

for n=1:number_of_elements

   vertices=P(:,T(:,n)); 
   vertices_zhongxin=sum(vertices')/3;  %%重心坐标
   K=Triangle_area(vertices);           %%根据三角形节点坐标计算面积
   T_tri=T(:,n);                        %%三角形网格的节点编号
%    for i=1:length(T_tri)
       
%    int_value=-feval(coe_fun,vertices_zhongxin(1),vertices_zhongxin(2))/4/K*(ai*aj+bi*bj);
   for alpha=1:number_of_local_basis_fun_trial
        for beta=1:number_of_local_basis_fun_test
            basis_index_trial=alpha
            basis_index_test=beta
            if alpha~=beta
            int_value=-K*Gauss_quad_2D_trial_test(coe_fun,1,vertices_zhongxin,vertices,basis_type_trial,basis_index_trial,basis_type_test,basis_index_test,[1 0],[1 0])+...
                      Gauss_quad_2D_trial_test(coe_fun,1,vertices_zhongxin,vertices,basis_type_trial,basis_index_trial,basis_type_test,basis_index_test,[0 1],[0 1]);

            
            j=Tb_test(beta,n);
            i=Tb_trial(alpha,n);
            A(i,j)=A(i,j)-int_value;
            A(i,i)=A(i,i)+int_value;
            end
        end
    end
end


end