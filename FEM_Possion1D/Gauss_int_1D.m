function [solution]=Gauss_int_1D(f,a,b,N)
%%%%f��һ���������
%%%%a�ǻ�������
%%%%b�ǻ�������
%%%%N�Ǹ�˹���ֵ����Ŀ
%%%%N��2��4��8��ȡֵ��ͨ��ȡ��8���ɴﵽ�㹻�ľ���
if N~=2&&N~=4&&N~=8
    warning='N�����������2��4��8��ȡֵ'
end
[aai,xxi]=generate_Gauss_local_1D([a,b],N);
solution=aai*f(xxi)';
end