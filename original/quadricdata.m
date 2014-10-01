function [x,y,z] = quadricdata(a,b,theta,d)

  % Create principal quadric...

  [lm_sx,lm_sy] = meshgrid(linspace(-0.5,0.5,d),linspace(-0.5,0.5,d));

  lm_sz = a.*lm_sx.^2 + b.*lm_sy.^2;

  % Flatten, remove all -negative z pts, cut out section -1<x<1

  ls_pts=0;
   for x=1:size(lm_sx,1)
     for y=1:size(lm_sy,2)
       %if (lm_sz(x,y)>0) 
         ls_pts=ls_pts+1;
         lm_nx(ls_pts) = lm_sx(x,y);
         lm_ny(ls_pts) = lm_sy(x,y);
         lm_nz(ls_pts) = lm_sz(x,y);
       %end;
     end;
   end;

   % Rotate points around z axis by theta (0->2pi)

   R = [cos(theta) -sin(theta) 0;...
       sin(theta)  cos(theta)  0;...
       0               0       1];

   pts = [lm_nx;lm_ny;lm_nz];
   pts = pts' * R;
   lm_nx = pts(:,1);
   lm_ny = pts(:,2);
   lm_nz = pts(:,3);

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

