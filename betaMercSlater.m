function [D] = betaMercSlater( Nr, Nc, numGr, obstGrid, epsilonR )
% fix epsilonC
if epsilonR == 0
  epsilonC = 1;
else
  epsilonC = 0;
end
% Find free sites
freeSites = find( obstGrid == 0 );
numFree = length( freeSites );
numObst = numGr - numFree;
% Build T matrix
Te = zeros( numFree, numFree );
Tl = zeros( numFree, numFree );
vl = zeros( numFree, 1 );
ve = zeros( numFree, 1 );
% loop over free site to build T matrix
for ii = 1:numFree
  % find indices
  [rT, cT] = ind2sub( [Nr, Nc], freeSites(ii) );
  % check to the right
  rIndObst = sub2ind( [Nr, Nc], rT, mod( cT + 1 - 1, Nc ) + 1);
  if obstGrid(rIndObst) == 0 %empty
    rIndT = find( freeSites == rIndObst );
    Tl(ii,rIndT) =  Tl(ii,rIndT) + 1 / 4;
    if epsilonC ~= 0
      Te(ii,rIndT) =  Te(ii,rIndT) - 1 / 4;
    end
  else
    Tl(ii,ii) = Tl(ii,ii) + 1 / 4;
    if epsilonC ~= 0
      Te(ii,ii) = Te(ii,ii) + 1 / 4;
    end
  end
  % check to the left
  lIndObst = sub2ind( [Nr, Nc], rT, mod( cT - 1 -1, Nc ) + 1);
  if obstGrid(lIndObst) == 0 %empty
    lIndT = find( freeSites == lIndObst );
    Tl(ii,lIndT) = Tl(ii,lIndT) + 1 / 4;
    if epsilonC ~= 0
      Te(ii,lIndT) = Te(ii,lIndT) + 1 / 4;
    end
  else
    Tl(ii,ii) = Tl(ii,ii) + 1 / 4;
    if epsilonC ~= 0
      Te(ii,ii) = Te(ii,ii) - 1 / 4;
    end
  end
  % check to below
  bIndObst = sub2ind( [Nr, Nc], mod( rT + 1 - 1, Nr ) + 1,  cT ) ;
  if obstGrid(bIndObst) == 0 %empty
    bIndT = find( freeSites == bIndObst );
    Tl(ii,bIndT) = Tl(ii,bIndT) + 1 / 4;
    if epsilonR ~= 0
      Te(ii,bIndT) = Te(ii,bIndT) - 1 / 4;
    end
  else
    Tl(ii,ii) = Tl(ii,ii) + 1 / 4;
    if epsilonR ~= 0
      Te(ii,ii) = Te(ii,ii) + 1 / 4;
    end
  end
  % check to above
  uIndObst = sub2ind( [Nr, Nc], mod( rT - 1 -1, Nr ) + 1,  cT );
  if obstGrid(uIndObst) == 0 %empty
    uIndT = find( freeSites == uIndObst );
    Tl(ii,uIndT) = Tl(ii,uIndT) + 1 / 4;
    if epsilonR ~= 0
      Te(ii,uIndT) = Te(ii,uIndT) + 1 / 4;
    end
  else
    Tl(ii,ii) = Tl(ii,ii) + 1 / 4;
    if epsilonR ~= 0
      Te(ii,ii) = Te(ii,ii) - 1 / 4;
    end
  end
  % velocity
  if epsilonR ~= 0
    vl(ii) = 1 / 4 * ( (1 - obstGrid(bIndObst) ) - (1 - obstGrid(uIndObst) ) );
    ve(ii) = 1 / 4 * ( (1 - obstGrid(bIndObst) ) + (1 - obstGrid(uIndObst) ) );
  else
    vl(ii) = 1 / 4 * ( (1 - obstGrid(rIndObst) ) - (1 - obstGrid(lIndObst) ) );
    ve(ii) = 1 / 4 * ( (1 - obstGrid(rIndObst) ) + (1 - obstGrid(lIndObst) ) );
  end
  %keyboard
end
% Solve matrix equation
Al = Tl-eye(numFree);
Al(end,:) = 1;
Ae = Te;
b = zeros(numFree,1);
b(end) = 1;
nl = linsolve( Al, b );
ne = linsolve( Al, -Ae * nl );
% Calculate D
D = ve' * nl + vl' * ne;
% Scale it
D = 2 * D;
keyboard
%% from the paper
% Tl_paper = ...
%   [ 1/2 1/4 0 1/4 0 0;...
%   1/4 1/2 1/4 0 0 0;...
%   0 1/4 0 1/4 1/2 0;...
%   1/4 0 1/4 0 0 1/2;...
%   0 0 1/2 0 1/4 1/4;...
%   0 0 0 1/2 1/4 1/4];
% Te_paper = ...
%   [ 0 1/4 0 -1/4 0 0;...
%   -1/4 0 1/4 0 0 0;...
%   0 -1/4 0 1/4 0 0;...
%   1/4 0 -1/4 0 0 0;...
%   0 0 0 0 1/4 1/4;...
%   0 0 0 0 -1/4 -1/4];
% 
% Te_paper2 = ...
%   [ 0 -1/4 0 1/4 0 0;...
%   1/4 0 -1/4 0 0 0;...
%   0 1/4 0 -1/4 0 0;...
%   -1/4 0 1/4 0 0 0;...
%   0 0 0 0 -1/4 -1/4;...
%   0 0 0 0 1/4 1/4];
% 
% vl_paper = [0; 0; 0; 0; -1/4; 1/4];
% vl_paper2 = [0; 0; 0; 0; 1/4; -1/4];
% ve_paper = [1/2; 1/2; 1/2; 1/2; 1/4; 1/4];
% 
% 
% Al_paper = Tl_paper-eye(numFree);
% Al_paper(end,:) = 1;
% Ae_paper = Te_paper;
% Ae_paper2 = Te_paper2;
% nl_paper = linsolve( Al_paper, b );
% ne_paper = linsolve( Al_paper, -Ae_paper * nl_paper );
% ne_paper2 = linsolve( Al_paper, -Ae_paper2 * nl_paper );
% 
% D_paper =  ve_paper' * nl_paper + vl_paper' * ne_paper;
% D_paper2 =  ve_paper' * nl_paper + vl_paper2' * ne_paper2;
% 
% nl_test = [ 1/6; 1/6; 1/6; 1/6; 1/6; 1/6];
% 
