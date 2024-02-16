function [y0,y1,y2,y3,y4,y5,y6] = bidwa_function(u0,u1,u2,u3,u4,u5,u6,rot_ind0,rot_ind1,direction)

if(direction == 0)
    U = double( [u0 u1 u2 u3 u4 u5 u6] );
    Y = circshift(U,[0 rot_ind0]);
else
    U = double( [u6 u5 u4 u3 u2 u1 u0] );
    Y = circshift(U,[0 -rot_ind1]);
end
y0 = Y(1);
y1 = Y(2);
y2 = Y(3);
y3 = Y(4);
y4 = Y(5);
y5 = Y(6);
y6 = Y(7);