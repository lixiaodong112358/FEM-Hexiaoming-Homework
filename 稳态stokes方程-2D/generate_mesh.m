function mesh=generate_mesh(N1,N2,basis_type)
%%%��������ڵ�ȫ�ֱ�ŵľ���2021/5/6��������
if basis_type==201
    KK=(N1+1)*(N2+1);
    RR=N2+1;
    CC=N1+1;
elseif basis_type==202
    KK=(2*N1+1)*(2*N2+1);
    RR=2*N2+1;
    CC=2*N1+1;
else
    warning='basis_type �������'
end
mesh=flipud(reshape([1:KK],RR,CC));
end