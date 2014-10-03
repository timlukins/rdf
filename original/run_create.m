function [X,Y,E,S,R] = run_create(type,varargin)

tt = 10;   % Number of frames 
q  = 1.0;
r  = 0.5;  % Radius (if used)
s  = 0.5;
d  = 101;  % Size of depth data
w  = 10;
mx = 0; % Translation movement
my = 0;
file = ''; % test.avi
depth = 0; % Generate tif depthmap output files
output = 1; % Generate tif depthmap/colour output files

if (length(varargin)>=1) 
	tt = varargin{1}; 
end;
if (length(varargin)>=2) 
	r = varargin{2}; 
end;
if (length(varargin)>=3) 
	d = varargin{3}; 
end;
if (length(varargin)>=4) 
	file = varargin{4}; 
end;

img = imread('brick_wall.jpg');
img = double(img)./255.0;

thetas = linspace(0,pi/6,tt);
radii = linspace(0.3,1.0,tt);
stretch = linspace(0.3,1.0,tt);
warp = linspace(-1.0,1.0,tt);
deltas = linspace(61,101,tt); % Stretch in x
gammas = linspace(61,101,tt); % Stretch in y
dispx = linspace(0,0.5,tt);
dispy = linspace(0,0,tt);

H=figure;
set(H,'Position',[100 100 600 600]);
set(H,'DoubleBuffer','on');
set(H,'Color','white');

if (~strcmp(file,''))
  disp(sprintf('*** WILL RECORD MOVIE TO %s',file));
  mov = avifile(file,'fps',1);
end;

lI = zeros(2);  % Last values
lII = zeros(2);
lIII = zeros(2);

X = zeros(tt,1);
Y = zeros(tt,1);
E = zeros(tt,1); % Values
S = zeros(tt,1);
R = zeros(tt,1);

for t=1:tt
   
%deltas(1)
%deltas(2)
% 
%    if (deltas(1)~=deltas(2))
%        deltas
%     px = (deltas(t)+1)/2; % Point on surface we are interested in
%     py = (gammas(t)+1)/2;
%    else
%     px = (deltas(1)+1)/2;
%     py = (gammas(1)+1)/2;
%    end;
px = (d+1)/2;
py = (d+1)/2;    

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  switch(lower(type))
    case 'rotcyl' 
      deltas = linspace(101,101,tt); 
      gammas = linspace(101,101,tt); 
      [x,y,z,f] = cylindata(r,thetas(t),d,mx,my);
    case 'movcyl'
      [x,y,z,f] = cylindata(r,0,d,dispx(t),my);  
    case 'expcyl'
      [x,y,z,f] = cylindata(radii(t),0,deltas(t),mx,my);
    case 'expsph'
      deltas = linspace(101,101,tt); 
      gammas = linspace(101,101,tt); 
      [x,y,z,f] = spheredata(radii(t),0,deltas(t));  
      deltas = linspace(101,61,tt); 
      gammas = linspace(101,61,tt);       
    case 'rexcyl'
      [x,y,z,f] = cylindata(radii(t),thetas(t),d);
    case 'roteli'
      [x,y,z,f] = ellipdata(q,r,s,thetas(t),d);
    case 'streli'
      deltas = linspace(101,61,tt); 
      gammas = linspace(101,101,tt); 
      [x,y,z,f] = ellipdata(stretch(t),r,s,0,d);      
    case 'rotqua' 
      [x,y,z,f] = quadricdata(q-1.5,r,thetas(t),d);
    case 'wrpqua' 
      [x,y,z,f] = quadricdata(q,warp(t),0,d);
    case 'random'
      [x,y] = meshgrid(linspace(-1,1,d),linspace(-1,1,d));
      %z = rand(d,d).*0.1;
      at = rand*6-3;
      [nx,ny]= meshgrid(linspace(at,at+0.3,d),linspace(at,at+0.3,d));
      z = peaks(nx,ny)*0.5; 
    otherwise
      disp('Usage:  run_create(type[,totaltime,radius,density,avifile]');
      disp('One of: rotcyl,expcyl,rexcyl,roteli,streli,rotqua,wrpqua, + random');
      disp('Where:  rot=rotate exp=expand rex=rotate&expand str=stretch wrp=warp');
      disp('And:    cyl=cylinder eli=elipsoid qua=quadric');
      close(gcf);
      return;
  end;
  
  %z = z + 0.001; % Not quite zero
  
  
  m = (z~=0); % Mask of moving points (simply > 0)
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$
  
  % Apply flow to image warp
  
  %f(:,2) = - f(:,2); % Image co-ords are -ve 
  %f = f .* size(z,1);
  %f = f + (size(z,1)/2); % Crucial rescale/translate to image coords
    %plot(f(:,1),f(:,2),'mx'); drawnow; pause(2); continue;  
  if (t==1)
      lastf = f; 
      %img = imresize(img,size(z)); % IF IMAGE ALWAYS BIGGER - DON' INCLUDE
      img = imresize(img,[size(z,1)*1.5,size(z,2)*1.5]); % 1.5 times bigger!
      % to accomodate for rotational/motion effects
      %img = img(1:size(z,1),1:size(z,2),:);
      %img = img(1:size(z,1)+101,1:size(z,2)+101,:);
  end;

  %quiver(lastf - f)
  T = cp2tform(lastf,f,'affine');%'piecewise linear');
  [wimg,xdata,ydata] = imtransform(img,T,'bicubic'); % Apply new transform to update original image

    %min(xdata)
    %min(ydata)
  %imshow(wimg); drawnow; pause(2); continue;  
  % with rotation/expansion need to cut out "centre" of warped image!
  cropfromx = floor((size(wimg,1)/2) - (size(z,1)/2));% + min(xdata);
  cropfromy = floor((size(wimg,2)/2) - (size(z,2)/2));% + min(ydata);
  wimg = wimg(cropfromx+1:cropfromx+size(z,1),cropfromy+1:cropfromy+size(z,2),:);
  %wimg = img;
    
  %aimg = m.*img;
  %bimg = m.*wimg;
  %wimg = aimg+bimg;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  dtx = deltas(t); % Sample grid size [-1,1]^2
  dty = gammas(t);
  
  % USING GRADIENT

%  [i,ii,iii] = fundforms(dtx,dty,x,y,z); 
%  I(:,:) = i(px,py,:,:); 
%  II(:,:) = ii(px,py,:,:); 
%  III(:,:) = iii(px,py,:,:); 
 
  % USING QUADRIC FITTING

  dx = x(px-w:py+w,px-w:py+w);
  dy = y(px-w:py+w,px-w:py+w); 
  dz = z(px-w:py+w,px-w:py+w); 
  nd = [reshape(dx,size(dx,1)*size(dx,2),1),reshape(dy,size(dy,1)*size(dy,2),1),reshape(dz,size(dz,1)*size(dz,2),1)];
    
  [I,II,III,nz] = fundforms(dtx,dty,nd); 
  
  [evec,eval]=eigs(II,2);
  k1 = abs(eval(1,1))*24; % SCALING CORRECTION!
  k2 = abs(eval(2,2))*24;
  r1 = evec(2,:); % OTHER WAY ROUND FOR SOME REASON
  r2 = evec(1,:);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  disp(sprintf('-------------------------------- step %d of %d',t,tt));

  %I
  %II
  %III

  if (t==1)
    lI = I;
    lII = II;
    D = [0 0; 0 0];
    P = [0 0; 0 0];
    lnr1 = r1;
    ldtx = dtx;
    ldty = dty;
  end;
  
  A = I-lI;
  B = II-lII;
  lI = I;
  lII = II;
  
  Q = D + inv([r1; r2])*B*P;
  D = [k1 0; 0 k2];
  P = [r1; r2];
  

%    [ls_r,ls_s,ls_e] = decomposechange(Q);
%    E(t) = ls_e;
%    S(t) = ls_s;
%    R(t) = ls_r;
  

  nr1 = r1;
  R(t) = acos(dot(lnr1,nr1))/pi;
  lnr1 = nr1;
  
  %e1 = B(1,1);
  %e2 = B(2,2);
  %R(t) = (abs(e1)+abs(e2))/2; 
  
  d1 = -(dtx - ldtx);%A(1,1)
  d2 = -(dty - ldty);%A(2,2)
  E(t) = (d1+d2)/2;
  ldtx = dtx;
  ldty = dty;
  
  if (d2==0)
    S(t) = d1;
  elseif (d1==0)
    S(t) = d2;
  elseif(d1>d2)
    S(t) = 1-(d2/d1);
  else
    S(t) = 1-(d1/d2);
  end;
  

  
% STORE EXTENT AND TYPE AS WELL!  
% 
%   case 'extent'
% 
%   % Work out EXTENT OF CHANGE (save only % inlier values within std) 
% 
%   lm_E  = sqrt(((lm_k1-lm_lk1).^2 + (lm_k2-lm_lk2).^2)/2);
%   lm_E = lm_E .* ((lm_E-mean(lm_E))<(std(lm_E)*ls_outreject));
%   %lm_sigmaE = lm_sigmaE/max(lm_sigmaE); 
%   %lm_sigmaE = log(lm_sigmaE);
% 
%   lm_sigma(f,:) = lm_sigma(f,:) + lm_E';
% 
% case 'deform'
% 
%   % Work out TYPE OF CHANGE
% 
%   lm_D = deformtype(lm_lk1,lm_lk2,(lm_lk1-lm_k1),(lm_lk2-lm_k2),25,25);
%   length(unique(lm_D))
% 
%   lm_sigma(f,:) = lm_sigma(f,:) + lm_D';%(lm_D~=lm_lD)';


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  clf;
  %figure;
  %surf(x,y,z,zeros(size(x)));
  %tri = delaunay(x,y);
  %trimesh(tri,x,y,z);
  surf(x,y,z,'Cdata',wimg,'EdgeColor','flat','FaceLighting','flat');
  hold on;
%  flow = f - lastf;
%size(flow)
%  quiver3(z,flow(:,1),flow(:,2),zeros(length(flow),1));
  %shading interp;
  %colormap(gray); 
  %surf(x,y,z);
  %hold on;
  plot3([0 r1(1)],[0 r1(2)],[z(px,py),z(px,py)],'r-');
  plot3([0 r2(1)],[0 r2(2)],[z(px,py),z(px,py)],'g-');
  %plot3(reshape(dx,size(dx,1)*size(dy,2),1),reshape(dy,size(dy,1)*size(dy,2),1),ones((w*2+1)^2,1)); 
  plot3(reshape(dx,size(dx,1)*size(dy,2),1),reshape(dy,size(dy,1)*size(dy,2),1),reshape(dz,size(dz,1)*size(dz,2),1),'r.'); 
  plot3(reshape(dx,size(dx,1)*size(dy,2),1),reshape(dy,size(dy,1)*size(dy,2),1),reshape(nz,size(nz,1)*size(nz,2),1),'b.'); 
  %plot(f(:,1),f(:,2),'mx');
  %legend('data','pdirect 1','pdirect 2','patch','quadric');  
  %colormap gray;
  axis equal;
  zlim([0 r*2]);
  view([45 30]); % OBLIQUE
  %view([0 90]); % TOP DOWN
  
  plot3([0 0],[0 0],[z(px,py),z(px,py)+1],'b-'); 
 

  title(sprintf('\\kappa_1 = %0.3f and \\kappa_2 = %0.3f',k1,k2));
  
  itext = sprintf(' %0.5f & %0.5f \\\\ %0.5f & %0.5f',I(1,1),I(1,2),I(2,1),I(2,2));
  %itext = strcat('$$ I = \left [ \begin{array}{cc} ',itext,' \end{array} \right ]$$');
  %itext = strcat('I = [ ',itext,' ] $$');
  %itext = sprintf('[E,F,G]=[%0.5f,%0.5f,%0.5f] \n [L,M,N]=[%0.5f,%0.5f,%0.5f] \n      A=[%0.5f,%0.5f,%0.5f] \n       B=[%0.5f,%0.5f,%0.5f]',I(1,1),I(1,2),I(2,2),II(1,1),II(1,2),II(2,2),A(1,1),A(1,2),A(2,2),B(1,1),B(1,2),B(2,2)); 
  itext = sprintf('[E,F,G]=[%0.5f,%0.5f,%0.5f]\n[L,M,N]=[%0.5f,%0.5f,%0.5f]',I(1,1),I(1,2),I(2,2),II(1,1),II(1,2),II(2,2)); 
  text('String', itext,...
       'Position',[1.0 -2.2 0],...
       'FontSize',12,...
       'FontNam','FixedWidth');
      %'Interpreter','latex',...
     
  drawnow;
  
  if (~strcmp(file,''))
    frame = getframe(gcf);
    mov = addframe(mov,frame);
  end;

  if (output==1)
    writepfm(double(z),sprintf('%03d.pfm',t)); 
    colour = uint8(round(wimg.*255));
    imwrite(colour,sprintf('%03d.ppm',t),'PPM');
    %writepfm(double(wimg),sprintf('colour_%03d.pbm',t));     
    %imwrite(double(z),sprintf('out_%d.pgm',t),'PGM'); 
    %imagesc(z);
  end;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end;

if (~strcmp(file,''))
   mov = close(mov);
end;

return;
