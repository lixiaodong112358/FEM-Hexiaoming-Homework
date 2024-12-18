function mesh=generate_mesh(N1,N2,basis_type)
%%%生成网格节点全局编号的矩阵，2021/5/6，李晓东
if basis_type==201
    KK=(N1+1)*(N2+1);
    RR=N2+1;
    CC=N1+1;
elseif basis_type==202
    KK=(2*N1+1)*(2*N2+1);
    RR=2*N2+1;
    CC=2*N1+1;
else
    warning='basis_type 输入错误'
end
mesh=flipud(reshape([1:KK],RR,CC));
end