function [X1,X0]=generate_initial_vector(Pb,dt)
%%����̬Stokes���̵ĳ�ʼֵ
%��ʼֵ
u1=feval('function_initial_u1',Pb(1,:),Pb(2,:));
u2=feval('function_initial_u2',Pb(1,:),Pb(2,:));
X0=[u1;u2];
%��һʱ��dt��λ��ֵ
%%����ʹ���˽��������ɳ�ֵ����ʵ���������Ҫ�ȽϺõĲ�ָ�ʽ����������Ϊ���˵��о���Ŀ��ʱ����Ҫ�񶯷��̣���ˣ���ʱ���Դ������о���
%%2021/09/11
u11=feval('fx1',Pb(1,:),Pb(2,:),dt,[0 0])';
u22=feval('fx2',Pb(1,:),Pb(2,:),dt,[0 0])';
% u11=u1+feval('function_initial_00_u1',Pb(1,:),Pb(2,:))'*dt;
% u22=u2+feval('function_initial_00_u2',Pb(1,:),Pb(2,:))'*dt;
X1=[u11;u22];
end