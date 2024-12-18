function result=local_FE_function_1D(x,uh_local_vec,NLB,vertices,basis_type,s)
%%%该函数测试通过，2021/4/30，李晓东；
result=0;
for k=1:NLB
    result=result+uh_local_vec(k)*FE_basis_local_fun_1D(x,vertices,basis_type,k,s);
end

end