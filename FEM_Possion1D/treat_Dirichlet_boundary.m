function [A,b]=treat_Dirichlet_boundary(function_b,Pb,A,b,boundarynodes)
%%%%����Dirichlet_boundary��0���������
nbn=size(boundarynodes,2);
%%
   for k=1:nbn
   if boundarynodes(1,k)==0
   
   i=boundarynodes(2,k);
   A(i,:)=0;
   A(i,i)=1;
   b(i)=feval(function_b,Pb(i));
   elseif boundarynodes(1,k)==1
       i=boundarynodes(2,k);
       b(i)=b(i)+exp(1)*(cos(1)-sin(1));
   end
   
end