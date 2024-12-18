 function result=reference_basis_2D(xi,yi,basis_type,basis_index,basis_der_x,basis_der_y)
 %201:2D linear
 %202:2D quad
 if basis_type==201
     if basis_index==1
         if basis_der_x==0&&basis_der_y==0
              result=-xi-yi+1;
          elseif basis_der_x==1&&basis_der_y==0
              result=-1;
          elseif basis_der_x==0&&basis_der_y==1
              result=-1;
         else
             warning='基函数的导数阶输入错误！'
         end
     elseif basis_index==2
         if basis_der_x==0&&basis_der_y==0
             result=xi;
         elseif basis_der_x==1&&basis_der_y==0
             result=1;
         elseif basis_der_x==0&&basis_der_y==1
             result=0;
         else
             warning='基函数的导数阶输入错误！'
         end
     elseif basis_index==3
         if basis_der_x==0&&basis_der_y==0
             result=yi;
         elseif basis_der_x==1&&basis_der_y==0
             result=0;
         elseif basis_der_x==0&&basis_der_y==1
             result=1;
         else
             warning='基函数的导数阶输入错误！'
         end
     else
         warning='basis_index输入错误,请在1-3中取值'
     end
 elseif basis_type==202
     if basis_index==1
          if basis_der_x==0&&basis_der_y==0
              result=2*xi^2+2*yi^2+4*xi*yi-3*yi-3*xi+1;
          elseif basis_der_x==1&&basis_der_y==0
              result=4*xi+4*yi-3;
          elseif basis_der_x==0&&basis_der_y==1
              result=4*xi+4*yi-3;
          elseif basis_der_x==2&&basis_der_y==0
              result=4;
          elseif basis_der_x==0&&basis_der_y==2
              result=4;
          elseif basis_der_x==1&&basis_der_y==1
              result=4;
          else
              warning='基函数的导数阶输入错误！'
          end
     elseif basis_index==2
         if basis_der_x==0&&basis_der_y==0
             result=2*xi^2-xi;
         elseif basis_der_x==1&&basis_der_y==0
             result=4*xi-1;
         elseif basis_der_x==0&&basis_der_y==1
             result=0;
         elseif basis_der_x==2&&basis_der_y==0
             result=4;
         elseif basis_der_x==0&&basis_der_y==2
             result=0;
         elseif basis_der_x==1&&basis_der_y==1
             result=0;
         else
             warning='基函数的导数阶输入错误！'
         end
     elseif basis_index==3
         if basis_der_x==0&&basis_der_y==0
             result=2*yi^2-yi;
         elseif basis_der_x==1&&basis_der_y==0
             result=0;
         elseif basis_der_x==0&&basis_der_y==1
             result=4*yi-1;
         elseif basis_der_x==2&&basis_der_y==0
             result=0;
         elseif basis_der_x==0&&basis_der_y==2
             result=4;
         elseif basis_der_x==1&&basis_der_y==1
             result=0;
         else
             warning='基函数的导数阶输入错误！'
         end
     elseif basis_index==4
         if basis_der_x==0&&basis_der_y==0
             result=-4*xi^2-4*xi*yi+4*xi;
         elseif basis_der_x==1&&basis_der_y==0
             result=4-8*xi-4*yi;
         elseif basis_der_x==0&&basis_der_y==1
             result=-4*xi;
         elseif basis_der_x==2&&basis_der_y==0
             result=-8;
         elseif basis_der_x==0&&basis_der_y==2
             result=0;
         elseif basis_der_x==1&&basis_der_y==1
             result=-4;
         else
             warning='基函数的导数阶输入错误！'
         end
     elseif basis_index==5
         if basis_der_x==0&&basis_der_y==0
             result=4*xi*yi;
         elseif basis_der_x==1&&basis_der_y==0
             result=4*yi;
         elseif basis_der_x==0&&basis_der_y==1
             result=4*xi;
         elseif basis_der_x==2&&basis_der_y==0
             result=0;
         elseif basis_der_x==0&&basis_der_y==2
             result=0;
         elseif basis_der_x==1&&basis_der_y==1
             result=4;
         else
             warning='基函数的导数阶输入错误！'
         end         
     elseif basis_index==6
         if basis_der_x==0&&basis_der_y==0
             result=-4*yi^2-4*xi*yi+4*yi;
         elseif basis_der_x==1&&basis_der_y==0
             result=-4*yi;
         elseif basis_der_x==0&&basis_der_y==1
             result=4-4*xi-8*yi;
         elseif basis_der_x==2&&basis_der_y==0
             result=0;
         elseif basis_der_x==0&&basis_der_y==2
             result=-8;
         elseif basis_der_x==1&&basis_der_y==1
             result=-4;
         else
             warning='基函数的导数阶输入错误！'
         end
     else
         warning='basis_index输入错误，请在1至6之间取值'
     end
 
 end