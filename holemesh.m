function [nels nnodes nodexy con bound left right bottom top void nvoidE nvoidN] = holemesh(L0, H0, f0, Nr, Nt)

% Input
% -----
% L0 - half-width of the rectangular domain
% H0 - half-height of the rectangular domain
% Nr - No of elements in the radial direction
% Nt - No of elements in the tangential direction
%
% Output
% ------
% nels - Total no of elements
% nnodes - Total no of nodes
% nodexy - Nodal coordinates
% con - Connectivity matrix
% bound - All nodes on the exterior boundaries
% left,right,bottom,top - Domain edges
% void - Void boundary segments


tol = 1.e-12;

R0 = sqrt(L0*H0);

r0 = sqrt(f0*4*L0*H0/pi);

bias = 1.1;

Ntn=Nt;
Nrn=Nr+1;
Nn=Ntn*Nrn;
Ne=Nt*Nr;
Ntr=Nrn*(Nt/8+1);
Ntl=Nrn*(3*Nt/8+1);
Nbl=Nrn*(5*Nt/8+1);
Nbr=Nrn*(7*Nt/8+1);

theta=2.*pi/Nt;
xcend=r0*cos(-theta);
ycend=r0*sin(-theta);
yrend=-2.*H0/(Nt/4);

x = [];
y = [];

n = 0;
for i=1:Nt
  x1 = r0*cos((i-1)*theta);
  y1 = r0*sin((i-1)*theta);

  if (x1 > tol)
    x2 = R0;
    y2 = (y1/x1)*x2;
    if (abs(y2) <= R0)
      xend = x2;
      yend = y2;
    end
  elseif (x1 < -tol)
    x2 = -R0;
    y2 = (y1/x1)*x2;
    if (abs(y2) <= R0)
      xend = x2;
      yend = y2;
    end
  else
    xend = 0;
    yend = R0*sign(y1);
  end

  if (y1 > tol)
    y2 = R0;
    x2 = (x1/y1)*y2;
    if (abs(x2) <= R0)
      xend = x2;
      yend = y2;
    end
  elseif (y1 < -tol)
    y2 = -R0;
    x2 = (x1/y1)*y2;
    if (abs(x2) <= R0)
      xend = x2;
      yend = y2;
    end
  else
    yend = 0;
    xend = R0*sign(x1);
  end

  if (bias ~= 1)
    xinc = (xend-x1)*(bias-1)/(bias^Nr-1);
    yinc = (yend-y1)*(bias-1)/(bias^Nr-1);
  else
    xinc = (xend-x1)/Nr;
    yinc = (yend-y1)/Nr;
  end

  x = [x x1];
  y = [y y1];
  for j=1:Nr
    if (bias ~= 1)
      x = [x x1+xinc*(bias^j-1)/(bias-1)*(L0-abs(x1))/(R0-abs(x1))];
      y = [y y1+yinc*(bias^j-1)/(bias-1)*(H0-abs(y1))/(R0-abs(y1))];
    else
      x = [x x1+j*xinc*(L0-abs(x1))/(R0-abs(x1))];
      y = [y y1+j*yinc*(H0-abs(y1))/(R0-abs(y1))];
    end
  end
end

el = [];

for i=1:Nt-1
for j=1:Nr
  n1 = (i-1)*Nrn+j;
  n2 = n1+1;
  n3 = n2+Nrn;
  n4 = n3-1;

  el = [el [n1;n2;n3;n4]];
end
end

for j=1:Nr
  n1 = (Nt-1)*Nrn+j;
  n2 = n1+1;
  n3 = j+1;
  n4 = j;

  el = [el [n1;n2;n3;n4]];
end

bound = [];
void = [];
right = [];
top = [];
left = [];
bottom = [];

for i=1:Nt-1
  void = [void [1+(i-1)*Nrn; 1+i*Nrn]];
end
void = [void [1+(Nt-1)*Nrn; 1]];

nvoidE = 1;
nvoidN = 1+Nt/4*Nrn;

n1 = Ntr;
n2 = Ntr+Nrn;
for i=1:Nt/4
  bound = [bound n1];
  top = [top [n1; n2]];
  n1 = n1+Nrn;
  n2 = n2+Nrn;
end

n1 = Ntl;
n2 = Ntl+Nrn;
for i=1:Nt/4
  bound = [bound n1];
  left = [left [n1; n2]];
  n1 = n1+Nrn;
  n2 = n2+Nrn;
end

n1 = Nbl;
n2 = Nbl+Nrn;
for i=1:Nt/4
  bound = [bound n1];
  bottom = [bottom [n1; n2]];
  n1 = n1+Nrn;
  n2 = n2+Nrn;
end

n1 = Nbr;
n2 = Nbr+Nrn;
for i=1:Nt/4
  bound = [bound n1];
  right = [right [n1; n2]];
  n1 = n1+Nrn;
  n2 = n2+Nrn;
  if (n1 > Nn)
    n1 = Nrn;
  end
  if (n2 > Nn)
    n2 = Nrn;
  end
end

nels = Ne;
nnodes = Nn;
nodexy = [x; y];
con = el;

% hold all
% axis equal
% 
% for i=1:length(el(1,:))
%   n1 = el(1,i);
%   n2 = el(2,i);
%   n3 = el(3,i);
%   n4 = el(4,i);
% 
%   plot([x(n1) x(n2)], [y(n1) y(n2)], 'b');
%   plot([x(n2) x(n3)], [y(n2) y(n3)], 'b');
%   plot([x(n3) x(n4)], [y(n3) y(n4)], 'b');
%   plot([x(n4) x(n1)], [y(n4) y(n1)], 'b');
% end
% 
% for i=1:length(right(1,:))
%   n1 = right(1,i);
%   n2 = right(2,i);
% 
%   plot([x(n1) x(n2)], [y(n1) y(n2)], 'r');
% end

end
