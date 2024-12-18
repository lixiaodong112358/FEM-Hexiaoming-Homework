function result=local_FE_function_2D(x,y,uh_local_vec,NLB,vertices,basis_type,s)
%%%该函数测试通过，2021/4/30，李晓东；
% NLB:number of local basis  基函数的数量
result=0;
for k=1:NLB
    result=result+uh_local_vec(k)*local_basis_2D(x,y,vertices,basis_type,k,s(1),s(2));
end

end