function [y0,y1,y2,y3,y4,y5,y6] = dwa_function(u0,u1,u2,u3,u4,u5,u6,rot_ind)
    

U = double( [u0 u1 u2 u3 u4 u5 u6] );

Y = circshift(U,[0 rot_ind]);

y0 = Y(1);
y1 = Y(2);
y2 = Y(3);
y3 = Y(4);
y4 = Y(5);
y5 = Y(6);
y6 = Y(7);