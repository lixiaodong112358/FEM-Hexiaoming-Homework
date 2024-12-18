function [Gauss_weight,Gauss_nodes]=generate_Gauss_local_1D(vertices,Gauss_type)
%%%积分上下限
lower_bound=vertices(1);
upper_bound=vertices(2);
%确定节点和系数
if Gauss_type==4
    Ai=[0.3478548451,0.3478548451,0.6521451549,0.6521451549];
    xi=[0.8611363116,-0.8611363116,0.3399810436,-0.3399810436];
elseif Gauss_type==8
    Ai=[0.1012285363,0.1012285363,0.2223810345,0.2223810345,0.3137066459,0.3137066459,0.3626837834,0.3626837834];
    xi=[0.9602898565,-0.9602898565,0.7966664774,-0.7966664774,0.5255324099,-0.5255324099,0.1834346425,-0.1834346425];
elseif Gauss_type==2
    Ai=[1,1];
    xi=[-1/sqrt(3),1/sqrt(3)];
else
    warning='高斯点数请在2、4、8中选取'
end
%计算实际节点和系数
Gauss_weight=(upper_bound-lower_bound)*Ai/2;
Gauss_nodes=(upper_bound-lower_bound)*xi/2+(upper_bound+lower_bound)/2;

end