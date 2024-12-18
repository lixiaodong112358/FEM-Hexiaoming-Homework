function [solution]=Gauss_int_1D(f,a,b,N)
%%%%f是一个函数句柄
%%%%a是积分下限
%%%%b是积分上限
%%%%N是高斯积分点的数目
%%%%N在2，4，8中取值，通常取到8即可达到足够的精度
if N~=2&&N~=4&&N~=8
    warning='N输入错误，请在2，4，8中取值'
end
[aai,xxi]=generate_Gauss_local_1D([a,b],N);
solution=aai*f(xxi)';
end