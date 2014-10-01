% FITQUADRIC - Fit quadric surface to 3D Z(X,Y) patch data. 
%
%       [quad,{error}] = dum_fitquadric(data,{options});
%
%       Pass data points in Z(x,y) square 2D matrix.
%       Evaluation performed over square domain (-1<=[x,y]<=1).
%       Returned vector is 10 co-efficients of implict quadric function. 
%
%       1      2     3       4      5      6      7     8     9    10
%       ax^2 + by^2 +cz^2 + 2hxy + 2gxz + 2fyz + 2ux + 2vy + 2wz + d
%
%       Can return MSE of distance between quadric and data.
%
%       Initial random solution is selected.
%
%       Options to refine fitting at cost to performance by (in order applied):
%
%         1.    'algebraic'  true [DEFAULT], false - Replace random with min Taubin solution.
%         2.    'genetic'    true, false [DEFAULT] - Minimise fit via GA. 
%         3.    'euclidean'  true, false [DEFAULT] - Minimise distance via non-linear simplex.
%
%       NOTE: If all of these are disabled, then only the random vector is returned. 
%
%       If data is a NON-SQUARE array of Mx3 then method treats it as data in R3.
%       In which case, the method will first align the data on a local co-ordinate frame.
%       Where the origin is on the centre point, and the normal from it is aligned to the Z-axis.
%
%       e.g.
%              
%               dx = 20; dy = 20;
%               [x,y] = meshgrid(linspace(-1,1,dx),linspace(-1,1,dy));
%               z = -2 * x.^2 + 2 * y.^2; % Elliptic paraboloid
%               q = dum_fitquadric(z,'algebraic',true,'euclidean',true);
%               figure; hold on;
%               surf(z);
%               surf(dum_quadric(dx,dy,q));
%       
% DUMLIB 1.0

%               
% Damn Useful Matlab Library (DUMLIB)
%
% Copyright (C) 2005 Tim Lukins 
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

function [rv_coeffs,varargout] = dum_fitquadric(pm_P,varargin)

  %if (sum(sum(pm_P<0.0))>0)
  %  disp('dum_fitquadric: patch data has negative values!');
  %end;

  V = pm_P;

  lc_type      = 'general';
  lf_show      = false;
  lf_align     = false;
  lf_algebraic = true;
  lf_euclidean = false;
  lf_genetic   = false;

  if (size(pm_P,1)~=3 && size(pm_P,2)==3) % If we have data in R3 - align it!
    lf_align = true;
    ls_n = sqrt(size(pm_P,1));
  end;

  if (mod(length(varargin),2)==0) % Process any extra args... 
    for a=1:2:length(varargin) 
      switch(lower(varargin{a}))
        case 'type'
          lc_type = varargin{a+1};
        case 'show'
          lf_show = varargin{a+1};
        case 'algebraic'
          lf_algebraic = varargin{a+1};
        case 'euclidean'
          lf_euclidean = varargin{a+1};
        case 'genetic'
          lf_genetic = varargin{a+1};
        otherwise 
          disp(sprintf('Unknown option: %s',varargin{a}));
      end; 
    end;
  end;
  
  Q = zeros(5,1);

  if (lf_align)
    
    % Fit plane to get initial normal...
    
    ls_mu = pm_P(((size(pm_P,1)+1)/2),:); % Centre of array - point to fit to
    %ls_mu = mean(pm_P,1);
    lm_A = [pm_P(:,1)-ls_mu(1) pm_P(:,2)-ls_mu(2) pm_P(:,3)-ls_mu(3)];
    [U,S,W] = svd(lm_A,0);
    [s, i] = min(diag(S));
    a = W(1,i);
    b = W(2,i);
    c = W(3,i);
    if (c<0) % Help initial estimate to point up!
      n = [-a -b -c];
    else
      n = [a b c];
    end;

    % Iteratively revise normal till close as possible...
  
    %while(abs(1-n(3))>0.0000000001)
    %while(abs(1-n(3))>0.0000000001 && n(3)<0)
    %while(n(3)<0)
    %lm_R = [linspace(-1,1,ls_n^2)'.^2 linspace(-1,1,ls_n^2)'.*linspace(-1,1,ls_n^2)' linspace(-1,1,ls_n^2)'.^2 linspace(-1,1,ls_n^2)' linspace(-1,1,ls_n^2)'];

    once = 1;
    while(n(3)~=1.0 || once)
%disp('.');
      n = n/norm(n);
      r3 = n;
      r1 = ((eye(3)-n*n')*[0 0 1]'); % Align to Z axis
      r1 = r1'/norm(r1);
      r2 = cross(r3,r1);
      R = [r1; r2; r3;];
      ls_mu = pm_P(((size(pm_P,1)+1)/2),:); % Centre of data patch is local centre
      pm_P = pm_P - (ones(size(pm_P,1),1)*ls_mu);
      pm_P = (R*pm_P')';
      V = pm_P; % Save rotated data for plot at the end

      % Test quadric ...

      lm_R = [pm_P(:,1).^2 pm_P(:,1).*pm_P(:,2) pm_P(:,2).^2 pm_P(:,1) pm_P(:,2)];
      lv_S = pm_P(:,3);
      %Q = inv(lm_R'*lm_R)*lm_R'*lv_S;
      %Q = pinv(lm_R)*lv_S;
      Q = lscov(lm_R,lv_S); 
      n = [-Q(4) -Q(5) 1]/(1 + Q(4)^2 + Q(5)^2);
      
      once = 0;
    end;
%disp('+');
   
    pm_P = reshape(pm_P(:,3),ls_n,ls_n); % Finally, turn into array

  end;
  
  ls_xs = size(pm_P,1);
  ls_ys = size(pm_P,2);
  lv_xx = linspace(-1,1,ls_xs); % Domain -1,1
  lv_yy = linspace(-1,1,ls_ys);

  if (strcmp(lc_type,'principal'))
    rv_coeffs = Q(1:3,:); % Last best fitting aligned quadric
    ls_err = minfit(rv_coeffs,pm_P); % Initial error
  
  else
    % OTHERWISE - of generalised type...

    rv_coeffs = rand(10,1); % Random to start (in case straight to GA fitting)...
    
    if (lf_algebraic)

      % Do intial guess with Taubin fitting...

      lm_H = zeros(10);
      lm_DH = zeros(10);
      for xi=1:ls_xs
        for yi=1:ls_ys
          x = lv_xx(xi);
          y = lv_yy(yi);
          z = pm_P(xi,yi);
          lv_h = [x^2,y^2,z^2,2*x*y,2*x*z,2*y*z,2*x,2*y,2*z,1]'; 
          lm_H = lm_H + lv_h*lv_h';   
          J = [[2*x     0      0  ];...
               [ 0     2*y     0  ];...
               [ 0      0     2*z ];...
               [2*y    2*x     0  ];...
               [2*z     0     2*x ];...
               [ 0     2*z    2*y ];...
               [ 2      0      0 ];...
               [ 0      2      0 ];...
               [ 0      0      2 ];...
               [ 0      0      0 ]];
          lm_DH = lm_DH + J*J';
        end;
      end;
     
      warning off;
      [lm_vects,lm_vals] = eigs(lm_H,lm_DH,10); % Solved by generalised eigenvector problem...
      %[lm_vects,lm_vals] = eig(full(lm_H),full(lm_DH),'chol'); % Solved by generalised eigenvector problem...
     
      rv_coeffs = lm_vects(1:10,size(lm_vals,1)); % Smallest eigenvalued solution...

    end;

    ls_err = minfit(rv_coeffs,pm_P); % Initial error - either Random or Taubin

    % Maybe do some search with a GA...

    if (lf_genetic)

      ls_p = 100; % Population 100
      ls_o = 10; % 10 cooeffs 
      ls_top = 0.1; % Elitism
      ls_mu = 0.05; % Mutation
      ls_maxg = 100; % max generations
      lm_pop = ones(ls_p,1)*rv_coeffs'; % Init to current solution
      lm_pop(round(ls_p/2)+1:ls_p,:) = lm_pop(round(ls_p/2)+1:ls_p,:) + rand(round(ls_p/2),10); 

      % Iterate per generation

      ls_best = 0.0;
      ls_deltabest = 0.0;
      lv_index = zeros(ls_p,1);
      lm_npop = zeros(ls_p,ls_o);
      for g=1:ls_maxg

        % Selection - rank in order and weight over total pop...

        for i=1:ls_p
          lv_besti(i) = minfit(lm_pop(i,:)',pm_P); % Distance found to target  
        end;
        [lv_besti,lv_index] = sort(lv_besti);
        ls_best = lv_besti(1);
        
        % Stop if last pop...

        if (g==ls_p) break; end;

        % Create next generation...

        for i=1:ls_p

          if (i<=round(ls_top*ls_p))
            
            % Elitism - copy over top n% - breed remainder
            
            lm_npop(i,:) = lm_pop(lv_index(i),:); 

          else

            % Selection - tournament based

            father = ceil(rand()*ls_p); 
            mother = ceil(rand()*ls_p); 
            
            %           - rank based
           
            %father = lv_index(round((((log(rand()+0.00001))+5)/5)*ls_p));
            %mother = lv_index(round((((log(rand()+0.00001))+5)/5)*ls_p));
            
            % Crossover - single point 

            lm_npop(i,:) = breed(lm_pop(father,:),lm_pop(mother,:));
                
            % Mutation - low prob of altering single value... 

            if (rand()<ls_mu)
              lm_npop(i,:) = mutate(lm_npop(i,:)); 
            end;

          end;
        
        end;

        lm_pop = lm_npop;

      end;

    rv_coeffs = lm_pop(lv_index(1),:);

    end;
  
  end; 
  
  ls_err = minfit(rv_coeffs,pm_P); % After any genetic search... 

  % Possible fine tune fitting with Euclidean minimization...

  if (lf_euclidean)
    
    lt_optoptions = optimset( ...
          'Display', 'off', ... % iter
          'FunValCheck', 'on',...
          'TolX', 0.00000000001, ...
          'TolFun', 0.00000000001, ...
          'MaxFunEvals', 1000, ...  % WAS 5000
          'MaxIter', 1000 ... % WAS 5000
          );

    [rv_coeffs,ls_err,ls_flag] = fminsearch(@minfit,rv_coeffs,lt_optoptions,pm_P);
  
  end;

  %X = dum_quadric(size(pm_P,1),size(pm_P,2),rv_coeffs); X((size(pm_P,1)+1)/2,(size(pm_P,2)+1)/2)
  
  varargout{1} = ls_err;

  if (lf_show)
    clf; plot3(V(:,1), V(:,2), V(:,3),'r+');
    W = dum_quadric(ls_n,ls_n,rv_coeffs);
    W = reshape(W,ls_n^2,1);
    hold on; plot3(V(:,1), V(:,2), W,'bx');
    legend('data','quadric');
    drawnow;
    %pause(0.5);
  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Objective/fitness function for Euclidean fit...

function rs_err = minfit(pv_coeffs,pm_P) % ps_err

  % Instantiate quadric from co-effs... 
  
  lm_X = dum_quadric(size(pm_P,1),size(pm_P,2),pv_coeffs);

  % If evaluation has moved from origin - return initial error! WE DEFINITELY DON'T WANT THIS

  %if (lm_X((size(pm_P,1)+1)/2,(size(pm_P,2)+1)/2)~=0)
  %  disp('Warning: tried to move origin!'); 
  %  rs_err = 99999999999999999999;%ps_err;
  %  return;
  %end;

  % Calculate distances between points...

  lv_dist = reshape(abs(pm_P-lm_X),size(pm_P,1)*size(pm_P,2),1); 

  %[lm_Y,lv_dist] = kdtree(pm_P,lm_X);

  % Error in distances... 

  rs_err = mean(lv_dist.^2);
 
  %rs_err = median(lv_dist);
  
  %rs_err = 0.0;
  %ls_c = std(lv_dist(1:round(0.5*size(lv_dist,1))));
  %for p=1:length(lv_dist)
  %  if (abs(lv_dist(p)) <= ls_c)
  %    rs_err = rm_err + ((ls_c^2)/6).*(1-(1-(lv_dist(p)/ls_c)^2)^3); % Tukey m-estimator
  %  else
  %    rs_err = rm_err + (ls_c^2)/6;
  %  end;
  %end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rm_new = breed(lm_father,lm_mother)

  ls_len = size(lm_father,2);
  ls_at = ceil(ls_len*rand()); % Point of crossover
  rm_new(1,:) = lm_father(1,:);
  rm_new(1,ls_at:ls_len) = lm_mother(1,ls_at:ls_len);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rm_new = mutate(lm_this)

  ls_len = size(lm_this,2);
  ls_at = ceil(ls_len*rand()); % Point of muatation 
  rm_new(1,:) = lm_this(1,:);
  rm_new(1,ls_at) = rm_new(1,ls_at)+rand();
 
  
