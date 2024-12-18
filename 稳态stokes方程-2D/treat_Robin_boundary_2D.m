function [w,R]=treat_Robin_boundary_2D(function_cq,function_cr,P,T,Pb,Tb,Ai_1D,Xi_1D,boundaryedges,number_of_local_basis_fun_trial,number_of_local_basis_fun_test,basis_type_trial,basis_type_test)
%-1:Dirichlet
%-2:Neuman
%-3:Robin
%-Ai_1D :一维高斯权重
%-Xi_1D :一维高斯节点
%% 预处理阶段
w=zeros(size(Pb,2),1);
R=sparse(size(Pb,2),size(Pb,2));
nbe=size(boundaryedges,2);

%% 
for k=1:nbe
    if boundaryedges(1,k)==-3
        nk=boundaryedges(2,k);
        vertices=P(:,T(:,nk));
        end_point=Pb(:,boundaryedges([3,4],k));
        [Ai,xi]=generate_Gauss_2D_line(end_point,Ai_1D,Xi_1D);
        %形成向量w
        for beta =1:number_of_local_basis_fun_test
             int_value=Gauss_int_line_test(function_cq,Ai,xi,vertices,basis_type_test,beta,0,0);
             w(Tb(beta,nk),1)=w(Tb(beta,nk),1)+int_value;            
        end
        
        
        %形成矩阵R
        for alpha=1:number_of_local_basis_fun_trial
            for beta=1:number_of_local_basis_fun_test
                int_value=Gauss_int_line_trial_test(function_cr,Ai,xi,vertices,basis_type_trial,basis_type_test,alpha,beta,0,0);
                R(Tb(beta,nk),Tb(alpha,nk))=R(Tb(beta,nk),Tb(alpha,nk))+int_value;
            end
        end
    end
end


end