function [I,II,III,varargout]=fundforms(dx,dy,varargin)

  if (length(varargin)==3)
    x = varargin{1};
    y = varargin{2};
    z = varargin{3};
    doquadric = false;
    I = zeros(size(z,1),size(z,2),2,2);
    II = zeros(size(z,1),size(z,2),2,2);
    III = zeros(size(z,1),size(z,2),2,2);
  elseif (length(varargin)==1) % one 3xn point array
    d = varargin{1};
    if (size(d,1)~=3 && size(d,2)~=3)
      disp('Data is not Nx3 array in fundforms');
      return;
    end;
    doquadric = true;
    I = zeros(2,2);
    II = zeros(2,2);
    III = zeros(2,2);
  else
    disp('Incorrect data input to dum_fundforms');
    return;
  end;

  %%%%%%%%%%%%%%%%%%%%%%% work out gradients...

  if (doquadric)

    q = dum_fitquadric(d,'type','principal','euclidean',true);
    
    if (nargout==4)
      trans = d(((size(d,1)+1)/2),3); % Centre of array - point to fit to
      nd = dum_quadric(sqrt(length(d)),sqrt(length(d)),q);
      nd = nd + (ones(size(nd,1),size(nd,2))*trans); % Translate back to fitted point
      varargout(1) = {nd}; 
    end;
    
    a = q(1);
    b = q(2);
    c = q(3);

    % First and second order partial derivatives of
    %         2              2
    % Q := a x  + b x y + c y 

    x = dx;  % X & Y need to the be the "new length" of x, y
    y = dy;
    Xu = [x, 0, 2*a*x + b*y];
    Xv = [0, y, b*x + 2*c*y];
    E           =   dot(Xu,Xu);
    F           =   dot(Xu,Xv);
    G           =   dot(Xv,Xv);
   
    Xuu = [0 0 2*a];
    Xuv = [0 0 b];
    Xvv = [0 0 2*c];
    n           =   [0 0 1]; % Aligned along Z axis!
    L           =   dot(Xuu,n);
    M           =   dot(Xuv,n);
    N           =   dot(Xvv,n);
        
    I  = [E,F;F,G];
    II = [L,M;M,N];

    H = (E*N+G*L-2*F*M)/(2*(E*G-F^2));
    K = (L*N-M^2)/(E*G-F^2);

    III = 2*H*II - K*I; 

  else

    [xu,xv]     =   gradient(x);
    [xuu,xuv]   =   gradient(xu);
    [xvu,xvv]   =   gradient(xv);

    [yu,yv]     =   gradient(y);
    [yuu,yuv]   =   gradient(yu);
    [yvu,yvv]   =   gradient(yv);

    [zu,zv]     =   gradient(z);
    [zuu,zuv]   =   gradient(zu);
    [zvu,zvv]   =   gradient(zv);

    for i=1:(size(z,1))
      for j=1:(size(z,2))
        Xu          =   [xu(i,j) yu(i,j) zu(i,j)];
        Xv          =   [xv(i,j) yv(i,j) zv(i,j)];
        Xuu         =   [xuu(i,j) yuu(i,j) zuu(i,j)];
        Xuv         =   [xuv(i,j) yuv(i,j) zuv(i,j)];
        Xvv         =   [xvv(i,j) yvv(i,j) zvv(i,j)];
        E           =   dot(Xu,Xu);
        F           =   dot(Xu,Xv);
        G           =   dot(Xv,Xv);
        m           =   cross(Xu,Xv);
        n           =   m/sqrt(sum(m.*m));
        L           =   dot(Xuu,n);
        M           =   dot(Xuv,n);
        N           =   dot(Xvv,n);

        I(i,j,:,:)  = [E,F;F,G];
        II(i,j,:,:) = [L,M;M,N];

        H = (E*N+G*L-2*F*M)/(2*(E*G-F^2));
        K = (L*N-M^2)/(E*G-F^2);

        III(i,j,:,:) = 2*H*II(i,j,:,:) - K*I(i,j,:,:); 

      end;
    end;
  end;
