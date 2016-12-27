% function betaMercSlater()
% Grid
Nc = 4;
Nr = 2;
numGr = Nr * Nc;
% Temporary epsilon
epsilonR = 0.1;
epsilonC = 0;
bHop = 0;
% Obstacle Mat
energyGrid = zeros( Nr, Nc );
energyGrid(2,1) = inf;
energyGrid(2,2) = inf;
obstGrid = zeros( Nr, Nc );
obstGrid( energyGrid ~= 0 ) = 1;
% Find free sites
freeSites = find( obstGrid == 0 );
numFree = length( freeSites );
numObst = numGr - numFree;
% Build T matrix
Te = zeros( numGr, numGr );
%Tl = zeros( numGr, numGr );
Tl = eye( numGr );
vl = zeros( numGr, 1 );
ve = zeros( numGr, 1 );
for ii = 1:numGr
  % find indices
  [rT, cT] = ind2sub( [Nr, Nc], ii );
  % Site occupation
  if obstGrid(ii) ~= 0 
    siteObst = 1;
  else
    siteObst = 0;
  end
  %% check to the right
  rIndObst = sub2ind( [Nr, Nc], rT, mod( cT + 1 - 1, Nc ) + 1);
  if siteObst == 1 && obstGrid(rIndObst) == 1 && bHop == 0
    probMoveTo = 0;
    probMoveAway = 0;
    probMoveStay = 1/4;
  else
    probMoveTo = 1/4 * min( exp( -( energyGrid(ii) - energyGrid(rIndObst) ) ), 1 );
    probMoveAway = 1/4 * min( exp( -( energyGrid(rIndObst)  - energyGrid(ii) ) ), 1 );
    probMoveStay = 1/4 * ...
      ( 1 - min( exp( -( energyGrid(rIndObst)  - energyGrid(ii) ) ), 1 ) );
  end
  Tl(ii,rIndObst) =  Tl(ii,rIndObst) + probMoveTo;
  Tl(ii,ii) = Tl(ii,ii) - probMoveAway;
  if epsilonC ~= 0 
    Te(ii,rIndObst) =  Te(ii,rIndObst) - probMoveTo;
    Te(ii,ii) = Te(ii,ii) +  probMoveStay;
  end
  %% check to the left
  lIndObst = sub2ind( [Nr, Nc], rT, mod( cT - 1 -1, Nc ) + 1);
  if siteObst == 1 && obstGrid(lIndObst) ~=0 && bHop == 0
    probMoveTo = 0;
    probMoveAway = 0;
    probMoveStay = 1/4;
  else
    probMoveTo = 1/4 * min( exp( -( energyGrid(ii) - energyGrid(lIndObst) ) ), 1 );
    probMoveAway = 1/4 * min( exp( -( energyGrid(lIndObst)  - energyGrid(ii) ) ), 1 );
    probMoveStay = 1/4 * ...
      ( 1 - min( exp( -( energyGrid(lIndObst)  - energyGrid(ii) ) ), 1 ) );
  end
  Tl(ii,lIndObst) =  Tl(ii,lIndObst) + probMoveTo;
  Tl(ii,ii) = Tl(ii,ii) - probMoveAway;
  if epsilonC ~= 0 
    Te(ii,lIndObst) =  Te(ii,lIndObst) + probMoveTo;
    Te(ii,ii) = Te(ii,ii) -  probMoveStay;
  end
  %% check to below
  bIndObst = sub2ind( [Nr, Nc], mod( rT + 1 - 1, Nr ) + 1,  cT ) ;
  if siteObst == 1 && obstGrid(bIndObst) == 1 && bHop == 0
    probMoveTo = 0;
    probMoveAway = 0;
    probMoveStay = 1/4;
  else
    probMoveTo = 1/4 * min( exp( -( energyGrid(ii) - energyGrid(bIndObst) ) ), 1 );
    probMoveAway = 1/4 * min( exp( -( energyGrid(bIndObst)  - energyGrid(ii) ) ), 1 );
    probMoveStay = 1/4 * ...
      ( 1 - min( exp( -( energyGrid(bIndObst)  - energyGrid(ii) ) ), 1 ) );
  end
  Tl(ii,bIndObst) =  Tl(ii,bIndObst) + probMoveTo;
  Tl(ii,ii) = Tl(ii,ii) - probMoveAway;
  if epsilonR ~= 0 
    Te(ii,bIndObst) =  Te(ii,bIndObst) - probMoveTo;
    %Te(ii,ii) = Te(ii,ii) +  probMoveAway;
    Te(ii,ii) = Te(ii,ii) +  probMoveStay;
  end
  %% check to above
  uIndObst = sub2ind( [Nr, Nc], mod( rT - 1 -1, Nr ) + 1,  cT );
  if siteObst == 1 && obstGrid(uIndObst) ~=0 && bHop == 0
    probMoveTo = 0;
    probMoveAway = 0;
    probMoveStay = 1/4;
  else
    probMoveTo = 1/4 * min( exp( -( energyGrid(ii) - energyGrid(uIndObst) ) ), 1 );
    probMoveAway = 1/4 * min( exp( -( energyGrid(uIndObst)  - energyGrid(ii) ) ), 1 );
    probMoveStay = 1/4 * ...
      ( 1 - min( exp( -( energyGrid(uIndObst)  - energyGrid(ii) ) ), 1 ) );
  end
  Tl(ii,uIndObst) =  Tl(ii,uIndObst) + probMoveTo;
  Tl(ii,ii) = Tl(ii,ii) - probMoveAway;
  if epsilonR ~= 0 
    Te(ii,uIndObst) =  Te(ii,uIndObst) + probMoveTo;
    %Te(ii,ii) = Te(ii,ii) -  probMoveAway;
    Te(ii,ii) = Te(ii,ii) -  probMoveStay;
  end
  % velocity
  if epsilonR ~= 0
    %vl(ii) = 1 / 4 * ( (1 - obstGrid(bIndObst) ) - );
    %ve(ii) = 1 / 4 * ( (1 - obstGrid(bIndObst) ) + (1 - obstGrid(uIndObst) ) );
    if siteObst == 1 && bHop == 0
      vl(ii) = 1 / 4 * ( ( 1 - obstGrid(bIndObst) ) - (1 - obstGrid(uIndObst) ) );
      ve(ii) = 1 / 4 * ( ( 1 - obstGrid(bIndObst) ) + (1 - obstGrid(uIndObst) ) );
    else
      vl(ii) = 1 / 4 * ( min( exp( -( energyGrid(bIndObst) - energyGrid(ii) ) ), 1 ) ...
      - min( exp( - ( energyGrid(uIndObst) - energyGrid(ii) ) ), 1 ) );
      ve(ii) = 1 / 4 * ( min( exp( -( energyGrid(bIndObst) - energyGrid(ii) ) ), 1 ) ...
      + min( exp( - ( energyGrid(uIndObst) - energyGrid(ii) ) ), 1 ) );
    end
  else
    if siteObst == 1
      vl(ii) = 1 / 4 * ( ( 1 - obstGrid(rIndObst) ) - (1 - obstGrid(lIndObst) ) );
      ve(ii) = 1 / 4 * ( ( 1 - obstGrid(rIndObst) ) + (1 - obstGrid(lIndObst) ) );
    else
      vl(ii) = 1 / 4 * ( min( exp( -( energyGrid(rIndObst) - energyGrid(ii) ) ), 1 ) ...
      - min( exp( - ( energyGrid(lIndObst) - energyGrid(ii) ) ), 1 ) );
      ve(ii) = 1 / 4 * ( min( exp( -( energyGrid(rIndObst) - energyGrid(ii) ) ), 1 ) ...
      + min( exp( - ( energyGrid(lIndObst) - energyGrid(ii) ) ), 1 ) );
    end
  end
end
% Solve matrix equation
Al = Tl-eye(numGr);
Al(end,:) = 1;
Ae = Te;
b = zeros(numGr,1);
b(end) = 1;
nl = linsolve( Al, b );
ne = linsolve( Al, -Ae * nl );
D = ve' * nl + vl' * ne;
disp(D)

%% from the paper

%% Full mat
Tl_a = ...
  [ 1/2 1/2 1/4 0 0 0 1/4 0;...
    0 1/4 0 0 0 0 0 0;...
    1/4 0 1/2 1/2 1/4 0 0 0;...
    0 0 0 1/4 0 0 0 0;...
    0 0 1/4 0 0 1/2 1/4 0;...
    0 0 0 1/4 1/2 1/4 0 1/4;
    1/4 0 0 0 1/4 0 0 1/2;...
    0 1/4 0 0 0 1/4 1/2 1/4];
  
 Te_a = ...
   [ 0 0 -1/4 0 0 0 1/4 0;...
    0 1/4 0 0 0 0 0 0;....
    1/4 0 0 0 -1/4 0 0 0;...
    0 0 0 -1/4 0 0 0 0;...
    0 0 1/4 0 0 0 -1/4 0;...
    0 0 0 0 0 -1/4 0 -1/4;...
    1/4 0 0 0 -1/4 0 0 0;...
    0 0 0 0 0 1/4 0 1/4 ];

Al_a = Tl_a-eye(8);
Al_a(end,:) = 1;
Ae_a = Te_a;
b = zeros(8,1);
b(end) = 1;

vl_a = [0; -1/4; 0; 1/4; 0; 1/4; 0; -1/4];
ve_a = [1/2; 1/4; 1/2; 1/4; 1/2; 1/4; 1/2; 1/4];
nl_a = linsolve( Al_a, b );
ne_a = linsolve( Al_a, -Ae_a * nl_a );
D_a = ve_a' * nl_a + vl_a' * ne_a;

disp(D_a)

