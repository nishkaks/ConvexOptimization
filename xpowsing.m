function varargout = xpowsing(x,mode)

global numf numg 

argout = 0;
f = 0;
if bitand(mode,1) 
  numf = numf + 1;
  argout = argout + 1;
  n = size(x,1); n4 = n/4;
  for i=1:n4
    z = x((i-1)*4+1:i*4);
    f = f + (z(1)+10*z(2))^2 + 5*(z(3)-z(4))^2 + (z(2)-2*z(3))^4 + ...
	10*(z(1)-z(4))^4;
  end
  varargout(argout) = {f};
end
if bitand(mode,2) 
  numg = numg + 1;
  argout = argout + 1;
  n = size(x,1); n4 = n/4;
  g = zeros(n,1);
  for i=1:n4
    z = x((i-1)*4+1:i*4);
    gpart = [2*(z(1)+10*z(2)) + 40*(z(1)-z(4))^3;...
      20*(z(1)+10*z(2)) + 4*(z(2)-2*z(3))^3; ...
      10*(z(3)-z(4)) - 8*(z(2)-2*z(3))^3; ...
      -10*(z(3)-z(4)) - 40*(z(1)-z(4))^3];
    g((i-1)*4+1:i*4) = gpart;
  end
  varargout(argout) = {g};
end

return;
