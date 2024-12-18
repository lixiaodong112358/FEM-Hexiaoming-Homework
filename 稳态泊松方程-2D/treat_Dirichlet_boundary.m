function [A,b]=treat_Dirichlet_boundary(function_b,Pb,A,b,boundarynodes)
%%%%����Dirichlet_boundary��-1���������
%%%%Newman������-2
%-1:Dirichlet
%-2:Neuman
%-3:Robin
nbn=size(boundarynodes,2);
%%
   for k=1:nbn
   if boundarynodes(1,k)==-1
   i=boundarynodes(2,k);
   A(i,:)=0;
   A(i,i)=1;
   b(i)=feval(function_b,Pb(1,i),Pb(2,i));
   elseif boundarynodes(1,k)==-2
       warning='Newman������û����,�������';
   end
end