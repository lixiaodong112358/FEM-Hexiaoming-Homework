function [Gauss_weights,Gauss_nodes]=generate_Gauss_2D_line(vertices,Gauss_coefficient_reference_1D,Gauss_point_reference_1D)

x1=vertices(1,1);
x2=vertices(1,2);
y1=vertices(2,1);
y2=vertices(2,2);

Gpn=length(Gauss_coefficient_reference_1D);
Gauss_nodes=zeros(2,Gpn);
Gauss_weights=zeros(1,Gpn);
if x1==x2
    %%%垂直直线
    a=min(y1,y2);  %%积分下限
    b=max(y1,y2);  %%积分上限
    [Gauss_coefficient_local_1D,Gauss_point_local_1D]=generate_Gauss_local_1D(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,a,b);
    Gauss_nodes(1,:)=zeros(1,Gpn)+x1;
    Gauss_nodes(2,:)=Gauss_point_local_1D;
    Gauss_weights=Gauss_coefficient_local_1D;
elseif y1==y2
    a=min(x1,x2);  %%积分下限
    b=max(x1,x2);  %%积分上限
    [Gauss_coefficient_local_1D,Gauss_point_local_1D]=generate_Gauss_local_1D(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,a,b);
    Gauss_nodes(1,:)=Gauss_point_local_1D;
    Gauss_nodes(2,:)=zeros(1,Gpn)+y1;
    Gauss_weights=Gauss_coefficient_local_1D;
else   
    a=min(x1,x2);  %%积分下限
    b=max(x1,x2);  %%积分上限
    [Gauss_coefficient_local_1D,Gauss_point_local_1D]=generate_Gauss_local_1D(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,a,b);
    Gauss_nodes(1,:)=Gauss_point_local_1D;
        slope=(y2-y1)/(x2-x1);
        J=sqrt(1+slope^2);
    Gauss_nodes(2,:)=y1+slope*(Gauss_point_local_1D-x1);
    
    Gauss_weights=Gauss_coefficient_local_1D*J;
    
end