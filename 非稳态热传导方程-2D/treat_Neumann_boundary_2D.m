function result=treat_Neumann_boundary_2D(Neumann_fun,t,P,T,Pb,Tb,Ai_1D,Xi_1D,boundaryedges,number_of_local_basis_fun_test,basis_type_test)
%-1:Dirichlet
%-2:Neuman
%-3:Robin
%-Ai_1D :一维高斯权重
%-Xi_1D :一维高斯节点
result=zeros(size(Pb,2),1);
nbe=size(boundaryedges,2);
for k=1:nbe
    if boundaryedges(1,k)==-2
        nk=boundaryedges(2,k);
        vertices=P(:,T(:,nk));
        end_point=Pb(:,boundaryedges([3,4],k));
        [Ai,xi]=generate_Gauss_2D_line(end_point,Ai_1D,Xi_1D);
        for beta =1:number_of_local_basis_fun_test
             int_value=Gauss_int_line_test(Neumann_fun,t,Ai,xi,vertices,basis_type_test,beta,0,0);
             result(Tb(beta,nk),1)=result(Tb(beta,nk),1)+int_value;            
        end
    end
end