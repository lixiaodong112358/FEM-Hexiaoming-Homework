function result=treat_Neumann_boundary_2D(Neumann_fun,P,T,Pb,Tb,boundaryedges,number_of_local_basis_fun_test,basis_type_test)
%-1:Dirichlet
%-2:Neuman
%-3:Robin

Gauss_type=4;
result=sparse(size(Pb,2),1);
nbe=size(boundaryedges,2);
[Gauss_coefficient_reference_1D,Gauss_point_reference_1D]=generate_Gauss_reference_1D(Gauss_type);
for k=1:nbe
    if boundaryedges(1,k)==-2
        nk=boundaryedges(2,k);
        vertices=P(:,T(:,nk));
        end_point=P(:,boundaryedges([3,4],k));
        [Ai,xi]=generate_Gauss_2D_line(end_point,Gauss_coefficient_reference_1D,Gauss_point_reference_1D);
        for beta =1:number_of_local_basis_fun_test
             int_value=Gauss_int_line_test(Neumann_fun,Ai,xi,vertices,basis_type_test,beta,0,0);
             result(Tb(beta,nk),1)=result(Tb(beta,nk),1)+int_value;            
        end
    end
end