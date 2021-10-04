% PHASE   A utility to generate boundary points for 2D phase plot.
% INPUTS: X,Y vectors of values for the x and y axes.
% OUTPUT: P, a 2xK matrix of initial values that can be passed to
%   an ODE solver (e.g., RK4).
function p=phase(x,y);
  miny=min(y);
  maxy=max(y);
  x=x(:)';
  y=y(:)';
  nx=length(x);
  ny=length(y);
  y=y(2:ny-1);
  ny=ny-2;
  p=[x;(miny+zeros(1,nx))];
  if ny>0 p=[p [(max(x)+zeros(1,ny));y]]; end
  p=[p [fliplr(x);(maxy+zeros(1,nx))]];
  if ny>0 p=[p [(min(x)+zeros(1,ny));fliplr(y)]]; end