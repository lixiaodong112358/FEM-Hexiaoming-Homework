function result=fx2(x,y,s)
%% λ��u2�Ľ�����
if s(1)==0&&s(2)==0
result=-2/3*x.*y.^3+2-pi*sin(pi*x);
elseif s(1)==1&&s(2)==0
    result=-2/3*y.^3-pi^2*cos(pi*x);
elseif s(1)==0&&s(2)==1
    result=-2*x.*y.^2;
else
    warning='������ĵ��������������ȡֵ0����1'
end

end