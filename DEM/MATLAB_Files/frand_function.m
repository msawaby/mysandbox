function [y0,y1,y2,y3,y4,y5,y6] = frand_function(u0,u1,u2,u3,u4,u5,u6)

U = double( [u0 u1 u2 u3 u4 u5 u6] );
Y = zeros(size(U));

% for i = 1:length(U),
%     while(1)
%         k = randi([1,7],1,1);
%         if U(k) ~= -1,
%             Y(i) = U(k);
%             U(k) = -1;
%             break;
%         end
%     end
% end

indexes = 1:7;
for i = 1:length(U),
     k = randsrc(1,1,indexes);
     Y(i) = U(k);
     indexes(indexes == k) = [];
end

y0 = Y(1);
y1 = Y(2);
y2 = Y(3);
y3 = Y(4);
y4 = Y(5);
y5 = Y(6);
y6 = Y(7);