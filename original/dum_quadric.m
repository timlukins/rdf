% QUADRIC - Evaluate quadric surface to 3D Z(x,y) data. 
%
%       quad = dum_quadric(x,y,coeffs);
%
%       Where coeffs vector is in the order:
%
%       1      2     3       4      5      6      7     8     9    10
%       ax^2 + by^2 +cz^2 + 2hxy + 2gxz + 2fyz + 2ux + 2vy + 2wz + d
%
%       Can also return the (number) TYPE of quadric: [quad,type] = dum_quadric(...)
%
%       (1)     Real ellipsoid             
%       (2)     Imaginary ellipsoid      
%       (3)     Hyperboloid of one sheet
%       (4)     Hyperboloid of two sheets   
%       (5)     Real quadric cone          
%       (6)     Imaginary quadric cone  
%       (7)     Elliptic paraboloid         
%       (8)     Hyperbolic paraboloid       
%       (9)     Real elliptic cylinder      
%       (10)    Imaginary elliptic cylinder  
%       (11)    Hyperbolic cylinder          
%       (12)    Real intersecting planes    
%       (13)    Imaginary intersecting planes
%       (14)    Parabolic cylinder           
%       (15)    Real parallel planes         
%       (16)    Imaginary parallel planes    
%       (17)    Coincident planes           
%
%       Numbers 1,3,4,7,8 are the five NON-DEGENERATE forms. 
% 
% DUMLIB 1.0

%               
% Damn Useful Matlab Library (DUMLIB)
%
% Copyright (C) 2005,2006 Tim Lukins 
% 
% This library is free software; you can redistribute it and/or
% modify it under the terms of the GNU Library General Public
% License as published by the Free Software Foundation; either
% version 2 of the License, or (at your option) any later version.
% 
% This library is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% Library General Public License for more details.
% 
% You should have received a copy of the GNU Library General Public
% License along with this library; if not, write to the
% Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
% Boston, MA  02110-1301, USA.

function [z,varargout] = dum_quadric(pm_x,pm_y,pv_coeffs)

  % Define domain...
  
  if (size(pm_x,1)==1 && size(pm_y,2)==1) % Scalar increments
    pv_x = linspace(-1,1,pm_x);
    pv_y = linspace(-1,1,pm_y);
  else
    disp('dum_quadric needs just the resolution!');
  end;

  % Check what type of quadric we have...
  
  if (length(pv_coeffs)==3) % Principal
      
    for xi=1:length(pv_x)
      for yi=1:length(pv_y)

        a = pv_coeffs(1);
        b = pv_coeffs(2); 
        c = pv_coeffs(3);

        x = pv_x(xi);
        y = pv_y(yi);

        z(xi,yi) = a*x^2 + b*x*y + c*y^2; % PRINCIPAL

      end;
    end;
    
  elseif (length(pv_coeffs)==10) % Generalised

    raytrace = false;

    % Instantiate generalised quadric from co-effs... 

    if (~raytrace)
      
      a = pv_coeffs(1); 
      b = pv_coeffs(2); 
      c = pv_coeffs(3);
      h = pv_coeffs(4);
      g = pv_coeffs(5);
      f = pv_coeffs(6);
      u = pv_coeffs(7);
      v = pv_coeffs(8);
      w = pv_coeffs(9);
      d = pv_coeffs(10);

      for xi=1:length(pv_x)
        for yi=1:length(pv_y)
          
          x = pv_x(xi);
          y = pv_y(yi);
      
          z1 = (-2*g*x - 2*w - 2*f*y + 2*(g^2*x^2 + 2*g*x*w + 2*g*x*f*y + w^2 + 2*w*f*y + f^2*y^2 - c*a*x^2 - c*b*y^2 - 2*c*u*x - 2*c*h*x*y - c*d - 2*c*v*y)^(1/2))/(2*c);
    
          z2 = (-2*g*x - 2*w - 2*f*y - 2*(g^2*x^2 + 2*g*x*w + 2*g*x*f*y + w^2 + 2*w*f*y + f^2*y^2 - c*a*x^2 - c*b*y^2 - 2*c*u*x - 2*c*h*x*y - c*d - 2*c*v*y)^(1/2))/(2*c);
    
          if (abs(z1)<abs(z2)) 
            z(xi,yi)=real(z1);
          else
            z(xi,yi)=real(z2);
          end;
                                                      
        end;
      end;
        
      %ls_xn = 2*(a.*x + h.*y + g.*z + u); % Full normals
      %ls_yn = 2*(h.*x + b.*y + f.*z + v);
      %ls_zn = 2*(g.*x + f.*y + c.*z + w);

    else % Ray trace...

      % NOTE: as if in matrix form (Glassner, Introduction to Ray-tracing, 1989, pp. 68-69)

      a = pv_coeffs(1); % a
      e = pv_coeffs(2); % e
      h = pv_coeffs(3); % h
      b = pv_coeffs(4); % b
      c = pv_coeffs(5); % c
      f = pv_coeffs(6); % f
      d = pv_coeffs(7); % d
      g = pv_coeffs(8); % g
      i = pv_coeffs(9); % i
      j = pv_coeffs(10); % j

      % Ray-trace orthogonally upward onto intersection with quadric patch...
     
      z = zeros(length(pv_x),length(pv_y));

      for x=1:length(pv_x)
        for y=1:length(pv_y)

          t = 1.0;
          
          xo = pv_x(x); % Origin of ray
          yo = pv_y(y);
          zo = 0;
          xd = 0;% Direction of ray (positive Z unit vector)
          yd = 0;
          zd = 1;

          Aq = a*xd^2 + 2*b*xd*yd + 2*c*xd*zd + e*yd^2 + 2*f*yd*zd + h*zd^2;
          Bq = 2*(a*xo*xd + b*(xo*yd+xd*yo) + c*(xo*zd+xd*zo) + d*xd + e*yo*yd + f*(yo*zd+yd*zo) + g*yd + h*zo*zd + i*zd); 
          Cq = a*xo^2 + 2*b*xo*yo + 2*c*xo*zo + 2*d*xo + e*yo^2 + 2*f*yo*zo + 2*g*yo + h*zo^2 + 2*i*zo + j;
          
          %Aq = h; % If just solving for z!
          %Bq = 2*(c*xo + f*yo + i); 
          %Cq = a*xo^2 + 2*b*xo*yo +  2*d*xo + e*yo^2 + 2*g*yo +  j;

          if (Aq==0)

            z(x,y) = -Cq/Bq;

          else

            if ((Bq.^2-4*Aq*Cq)<0) % Discriminant says: No intersection

              z(x,y) = 0.0;

            else
           
              t = (-Bq - sqrt(Bq^2-4*Aq*Cq))/(2*Aq); % TEST THIS ROOT FIRST!
              
              if (t<0.0)
                t = (-Bq + sqrt(Bq^2-4*Aq*Cq))/(2*Aq);
              end;
          
            end;

          end;
          
          z(x,y) = t; % Since orthogonal - just distance up
          %z(x,y) = lv_Ro + lv_Rd .* t; % Ray intersection with quadric at distance t

        end;
      end;
    end;

    % If requesting type as well...

    if (nargout==2)
      
      % See: http://www.geom.uiuc.edu/docs/reference/CRC-formulas/node61.html

      a = pv_coeffs(1);
      b = pv_coeffs(2);
      c = pv_coeffs(3);
      h = pv_coeffs(4);
      g = pv_coeffs(5);
      f = pv_coeffs(6);
      p = pv_coeffs(7);
      q = pv_coeffs(8);
      r = pv_coeffs(9);
      d = pv_coeffs(10);

      lm_e = [a h g; h b f; g f c];
      lm_E = [a h g p; h b f q; g f c r; p q r d];

      p3 = rank(lm_e);
      p4 = rank(lm_E);
      delta = det(lm_E);
      ksame = ksigns(eig(lm_e));
      Ksame = ksigns(eig(lm_E));
     
      type = 0;

      if     (p3==3 & p4==4 & delta<0 & ksame) type = 1;
      elseif (p3==3 & p4==4 & delta>0 & ksame) type = 2;
      elseif (p3==3 & p4==4 & delta>0 & ~ksame) type = 3;
      elseif (p3==3 & p4==4 & delta<0 & ~ksame) type = 4;
      elseif (p3==3 & p4==3 & ~ksame) type = 5;
      elseif (p3==3 & p4==3 & ksame)  type = 6;
      elseif (p3==2 & p4==4 & delta<0 & ksame) type = 7;
      elseif (p3==2 & p4==4 & delta>0 & ~ksame)  type = 8;
      elseif (p3==2 & p4==3 & ksame & ~Ksame) type = 9;
      elseif (p3==2 & p4==3 & ksame & Ksame)  type = 10;
      elseif (p3==2 & p4==3 & ~ksame)  type = 11;
      elseif (p3==2 & p4==2 & ~ksame)  type = 12;
      elseif (p3==2 & p4==2 & ksame)  type = 13;
      elseif (p3==1 & p4==3)  type = 14;
      elseif (p3==1 & p4==2 & ~Ksame) type = 15;
      elseif (p3==1 & p4==2 & Ksame) type = 16;
      elseif (p3==1 & p4==1)  type = 17;
      end;

      varargout(1) = {type}; 
    end;

  else
    disp('Unrecognised number of quadric coefficients!');
  end;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function same=ksigns(vals)

  pos = 0;
  neg = 0;
  zer = 0;
  same = false;

  for v=1:length(vals)
    if (vals(v)>0) pos=pos+1; end;
    if (vals(v)<0) neg=neg+1; end;
    if (vals(v)==0) zer=zer+1; end;   
  end;

  if (pos>0 & neg==0) same=true; end;
  if (neg>0 & pos==0) same=true; end;

