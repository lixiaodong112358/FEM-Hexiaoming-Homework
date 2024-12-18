function BN=generate_boundarynodes(N1,N2,basis_type,BJ)
%%%���ɱ߽�ڵ�����
%%%2021/5/6 ������
    if basis_type==202
        RR=2*N2+1;
        CC=2*N1+1;
    elseif basis_type==201
        RR=N2+1;
        CC=N1+1;
    end
    k1=[1:RR:RR*CC-RR+1];
    k2=[RR*CC-RR+2:RR*CC];
    k3=RR*CC-RR:-RR:RR;
    k4=[RR-1:-1:2];
    BN=[k1,k2,k3,k4];
%%%%����ı߽�������������Ҫ������Ӧ�ص���
    if BJ=='D'
        BN=[zeros(1,size(BN,2))-1;BN];
    end
    BN(1,2:N1)=BN(1,2:N1)-1;   %%�ײ�Neumann�߽�  ע������Dirichilet  �߽��������ȱ�֤
end