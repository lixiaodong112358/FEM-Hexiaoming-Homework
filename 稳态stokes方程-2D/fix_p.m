function [AA,bb]=fix_p(fxp,A,b,P,Nb)
   %%将第一个节点
for i=3:4
   A(2*Nb+i,:)=0;
   A(2*Nb+i,2*Nb+i)=1;
   b(2*Nb+i)=feval(fxp,P(1,i),P(2,i),[0 0]);
   AA=A;
   bb=b;
end
end