function run_create(type,varargin)

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

if (length(varargin)>=1) tt = varargin{1}; end;
if (length(varargin)>=2) r = varargin{2}; end;
if (length(varargin)>=3) d = varargin{3}; end;
if (length(varargin)>=4) file = varargin{4}; end;

img = imread('brick_wall.jpg');
img = double(img)./255.0;
%img = img(1:d,1:d,:); % Just cut what we need

thetas = linspace(0,pi/8,tt);
radii = linspace(0.1,1.0,tt);
stretch = linspace(0.1,1.0,tt);
warp = linspace(-1.0,1.0,tt);
deltas = linspace(101,101,tt); % Stretch in x
gammas = linspace(101,101,tt); % Stretch in y
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

for t=1:tt
   
deltas(1)
deltas(2)

   if (deltas(1)~=deltas(2))
    px = (deltas(t)+1)/2; % Point on surface we are interested in
    py = (gammas(t)+1)/2;
   else
    px = (deltas(1)+1)/2;
    py = (gammas(1)+1)/2;
   end;
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  switch(lower(type))
    case 'rotcyl' 
      [x,y,z,f] = cylindata(r,thetas(t),d,mx,my);
    case 'movcyl'
      [x,y,z,f] = cylindata(r,0,d,dispx(t),my);  
    case 'expcyl'
      [x,y,z,f] = cylindata(radii(t),0,deltas(t),mx,my);
    case 'expsph'
      [x,y,z,f] = spheredata(radii(t),0,d);      
    case 'rexcyl'
      [x,y,z,f] = cylindata(radii(t),thetas(t),d);
    case 'roteli'
      [x,y,z,f] = ellipdata(q,r,s,thetas(t),d);
    case 'streli'
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

min(xdata)
min(ydata)
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
 
  % USING GRADIENT

  %[i,ii,iii] = fundforms(x,y,z); 
  %I(:,:) = i(px,py,:,:); 
  %II(:,:) = ii(px,py,:,:); 
  %III(:,:) = iii(px,py,:,:); 
 
  % USING QUADRIC FITTING

  dx = x(px-w:py+w,px-w:py+w);
  dy = y(px-w:py+w,px-w:py+w); 
  dz = z(px-w:py+w,px-w:py+w); 
  nd = [reshape(dx,size(dx,1)*size(dx,2),1),reshape(dy,size(dy,1)*size(dy,2),1),reshape(dz,size(dz,1)*size(dz,2),1)];
  
  dtx = 2.0/deltas(t); % Sample grid size!
  dty = 2.0/gammas(t);
  [I,II,III,nz] = fundforms(dtx,dty,nd); 
  
  [evec,eval]=eigs(II,2);
  k1 = eval(1,1);
  k2 = eval(2,2);
  r1 = evec(2,:); % OTHER WAY ROUND FOR SOME REASON
  r2 = evec(1,:);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  disp(sprintf('-------------------------------- step %d of %d',t,tt));

  I
  II
  %III

  A = I-lI;
  B = II-lII;

  %[U,S,V] = svd(II);
  %U
  %S
  %V

  lI = I;
  lII = II;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  clf;
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
  legend('data','pdirect 1','pdirect 2','patch','quadric');  
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
  itext = sprintf('[E,F,G]=[%0.5f,%0.5f,%0.5f] \t [L,M,N]=[%0.5f,%0.5f,%0.5f]\n      A=[%0.5f,%0.5f,%0.5f] \t       B=[%0.5f,%0.5f,%0.5f]',I(1,1),I(1,2),I(2,2),II(1,1),II(1,2),II(2,2),A(1,1),A(1,2),A(2,2),B(1,1),B(1,2),B(2,2)); 
  text('String', itext,...
       'Position',[0.5 -2.5 0],...
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
