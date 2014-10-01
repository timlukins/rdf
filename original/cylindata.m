function [x,y,z,f] = cylindata(r,theta,d,mx,my)

  % Create side on cylinder...

  [lm_sz,lm_sy,lm_sx]=cylinder(r,d);

  lm_sx = (lm_sx*2*r)-r;

  % Flatten, remove all -negative z pts, cut out section -1<x<1

  ls_pts=0;
   for x=1:size(lm_sx,1)
     for y=1:size(lm_sy,2)
       if (lm_sz(x,y)>0) 
         ls_pts=ls_pts+1;
         lm_nx(ls_pts) = lm_sx(x,y);
         lm_ny(ls_pts) = lm_sy(x,y);
         lm_nz(ls_pts) = lm_sz(x,y);
       end;
     end;
   end;

   % Rotate/translate points around z axis by theta (0->2pi)

   R = [cos(theta) -sin(theta) 0;...
       sin(theta)  cos(theta)  0;...
       0               0       1];
   
   T = [mx my 0]';
   M = [R T; 0 0 0 1];
   pts = [lm_nx;lm_ny;lm_nz;ones(1,ls_pts)];  
   npts = (M * pts)';

   lm_nx = npts(:,1);
   lm_ny = npts(:,2);
   lm_nz = npts(:,3);
   
   % Build 2D flow description
   
   f = [lm_nx, lm_ny];

   % Resample orthogonally 
 
   warning off MATLAB:griddata:DuplicateDataPoints
   %[x,y]=meshgrid(linspace(-r,r,d),linspace(-r,r,d));
   [x,y]=meshgrid(linspace(-1,1,d),linspace(-1,1,d));
   z = griddata(lm_nx,lm_ny,lm_nz,x,y);
   
   % Remove nans
 
   for xi=1:size(z,1)
     for yi=1:size(z,2)
       if (isnan(z(xi,yi)))
         z(xi,yi)=0;
       end;
     end;
   end;

