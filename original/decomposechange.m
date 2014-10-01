function [E,S,R] = decomposechange(dS)

  A = [1 0 0 1; 1 0 0 -1; 0 1 1 0; 0 -1 1 0]; % A mapping matrix
  %A = [1 0 0 1; 1 0 0 -1; 0 -1 1 0; 0 1 1 0]; % A mapping matrix

  B = reshape(dS,4,1); % B result reshaped

  %%%%%%%%%%

  x = lscov(A',B); % Solve by Least squares %s = v/q' not so accurate

  x(4) = -x(4); % Last value always needs correction % WHY?

  %%%%%%%%%%

  E = abs( x(1) );  % Expansion

  %S = abs( x(2) - x(3) );
  S = abs( x(2) + sqrt(-1)*x(3) ); % Shear - note use imaginary number
  
  R = abs( x(4) ); % Rotation
  
