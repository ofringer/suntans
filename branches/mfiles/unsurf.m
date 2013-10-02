% 
% UNSURF
%  H = UNSURF(TRI,X,Y,Z) fills the triangles in TRI which contains
%  indices to the vertices in X and Y with the color in Z.
%
function h = unsurf(varargin)
  
  
  t=varargin{1};
  x=varargin{2};
  y=varargin{3};
  z=varargin{4};

  held=ishold;

  xf = x([1+t(:,1:3),1+t(:,1)]);
  yf = y([1+t(:,1:3),1+t(:,1)]);  
  
  if(nargin>4)
    patch(xf',yf',z',varargin{5:end});
  else
    patch(xf',yf',z');
  end
  
  if(~held)
    hold off;
  end