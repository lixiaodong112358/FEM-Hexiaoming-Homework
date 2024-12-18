function [X1,X0]=generate_initial_vector(Pb,dt)
%%非稳态Stokes方程的初始值
%初始值
u1=feval('function_initial_u1',Pb(1,:),Pb(2,:));
u2=feval('function_initial_u2',Pb(1,:),Pb(2,:));
X0=[u1;u2];
%第一时刻dt的位移值
%%这里使用了解析解生成初值，真实情况可能需要比较好的差分格式处理，但是因为本人的研究项目暂时不需要振动方程，因此，暂时不对此深入研究；
%%2021/09/11
u11=feval('fx1',Pb(1,:),Pb(2,:),dt,[0 0])';
u22=feval('fx2',Pb(1,:),Pb(2,:),dt,[0 0])';
% u11=u1+feval('function_initial_00_u1',Pb(1,:),Pb(2,:))'*dt;
% u22=u2+feval('function_initial_00_u2',Pb(1,:),Pb(2,:))'*dt;
X1=[u11;u22];
end