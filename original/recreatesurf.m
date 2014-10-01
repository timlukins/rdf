function z=recreatesurf(x,y,k1,k2,r1,r2)

  if (dot(r1,r2)~=0)
    disp('Principle directions not orthogonal!');
    return;
  end;

  H = (1/2)*(k1+k2);
  K = k1*k2;
  
  alpha = (H-sqrt(H^2+K))/2;
  gamma = (H+sqrt(H^2+K))/2;
  theta = acos(dot(r1,[1 0]));  % Angle between first principal and X axis!
  
  a = alpha*cos(theta)^2 + gamma*sin(theta)^2;
  b = (alpha-gamma)*cos(theta)*sin(theta);
  c = alpha*sin(theta)^2 + gamma*cos(theta)^2;
  
  z = a*x.^2 + b*x*y * c*y.^2;
