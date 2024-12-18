function [result]=generate_initial_vector(function_initial_u1,function_initial_u2,function_initial_p,Pb,P)
%%����̬Stokes���̵ĳ�ʼֵ

u1=feval(function_initial_u1,Pb(1,:),Pb(2,:))';
u2=feval(function_initial_u2,Pb(1,:),Pb(2,:))';
p=feval(function_initial_p,P(1,:),P(2,:))';
result=[u1;u2;p];
end