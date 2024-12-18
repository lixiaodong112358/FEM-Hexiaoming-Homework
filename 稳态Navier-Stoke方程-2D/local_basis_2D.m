function result=local_basis_2D(x,y,vertices,basis_type,basis_index,basis_der_x,basis_der_y)

xn1=vertices(1,1);
xn2=vertices(1,2);
xn3=vertices(1,3);
yn1=vertices(2,1);
yn2=vertices(2,2);
yn3=vertices(2,3);


J=(xn2-xn1)*(yn3-yn1)-(xn3-xn1)*(yn2-yn1);
xi=((yn3-yn1)*(x-xn1)-(xn3-xn1)*(y-yn1))/J;
yi=(-(yn2-yn1)*(x-xn1)+(xn2-xn1)*(y-yn1))/J;
if basis_type==202
    if basis_der_x==0&&basis_der_y==0
        result=reference_basis_2D(xi,yi,basis_type,basis_index,0,0);
    elseif basis_der_x==1&&basis_der_y==0
        result=(yn3-yn1)/J*reference_basis_2D(xi,yi,basis_type,basis_index,1,0)+...
               (yn1-yn2)/J*reference_basis_2D(xi,yi,basis_type,basis_index,0,1);
    elseif basis_der_x==0&&basis_der_y==1
        result=(xn1-xn3)/J*reference_basis_2D(xi,yi,basis_type,basis_index,1,0)+...
               (xn2-xn1)/J*reference_basis_2D(xi,yi,basis_type,basis_index,0,1);    
    elseif basis_der_x==2&&basis_der_y==0
        result=(yn3-yn1)^2/J^2*reference_basis_2D(xi,yi,basis_type,basis_index,2,0)+...
               (yn1-yn2)^2/J^2*reference_basis_2D(xi,yi,basis_type,basis_index,0,2)+...
               2*(yn3-yn1)*(yn1-yn2)/J/J*reference_basis_2D(xi,yi,basis_type,basis_index,1,1);   
    elseif basis_der_x==0&&basis_der_y==2
        result=(xn1-xn3)^2/J^2*reference_basis_2D(xi,yi,basis_type,basis_index,2,0)+...
               (xn2-xn1)^2/J^2*reference_basis_2D(xi,yi,basis_type,basis_index,0,2)+...
               2*(xn1-xn3)*(xn2-xn1)/J/J*reference_basis_2D(xi,yi,basis_type,basis_index,1,1);
    elseif basis_der_x==1&&basis_der_y==1
        result=(xn1-xn3)*(yn3-yn1)/J^2*reference_basis_2D(xi,yi,basis_type,basis_index,2,0)+...
               (xn2-xn1)*(yn1-yn2)/J^2*reference_basis_2D(xi,yi,basis_type,basis_index,0,2)+...
               (xn1-xn3)*(yn1-yn2)/J/J*reference_basis_2D(xi,yi,basis_type,basis_index,1,1)+...
               (xn2-xn1)*(yn3-yn1)/J/J*reference_basis_2D(xi,yi,basis_type,basis_index,1,1);
    else
        warning='基函数导数阶输入错误！'
    end
elseif basis_type==201
    if basis_der_x==0&&basis_der_y==0
        result=reference_basis_2D(xi,yi,basis_type,basis_index,0,0);
    elseif basis_der_x==1&&basis_der_y==0
        result=(yn3-yn1)/J*reference_basis_2D(xi,yi,basis_type,basis_index,1,0)+...
               (yn1-yn2)/J*reference_basis_2D(xi,yi,basis_type,basis_index,0,1);
    elseif basis_der_x==0&&basis_der_y==1
        result=(xn1-xn3)/J*reference_basis_2D(xi,yi,basis_type,basis_index,1,0)+...
               (xn2-xn1)/J*reference_basis_2D(xi,yi,basis_type,basis_index,0,1);
    else
        warning='基函数导数阶输入错误！'
    end
end