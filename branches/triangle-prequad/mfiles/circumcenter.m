% CIRCUMCENTER Compute the circumcenter of a triangle.
%   [XV,YV]=CIRCUMCENTER(X1,Y1,X2,Y2,X3,Y3) computes the circumcenter
%   of the triangle defined by the vertices (x1,y1), (x2,y2), (x3,y3).

% Sheng Chen
% ICME, Stanford University
% Febrary 02, 2006
function [xv, yv] = circumcenter(x1, y1, x2, y2, x3, y3);

  xv = (x2.^2.*y1 - x3.^2.*y1 - x1.^2.*y2 + x3.^2.*y2 - y1.^2.*y2 + y1.*y2.^2 + x1.^2.*y3 ...
        - x2.^2.*y3 + y1.^2.*y3 - y2.^2.*y3 - y1.*y3.^2 + y2.*y3.^2)./ ...
       (2*(x2.*y1 - x3.*y1 - x1.*y2 + x3.*y2 + x1.*y3 - x2.*y3));
  yv = -(-(-2*x2 + 2*x3).*(x1.^2 - x2.^2 + y1.^2 - y2.^2) + ...
         (-2*x1 + 2*x2).*(x2.^2 - x3.^2 + y2.^2 - y3.^2))./ ...
       (-4*x2.*y1 + 4*x3.*y1 + 4*x1.*y2 - 4*x3.*y2 - 4*x1.*y3 + 4*x2.*y3);

