% Grid
Nc = 10; % number columns
Nr = 10; % number rows
numRuns = 10; % number of runs
% Temporary epsilon
epsilonR = 1;
epsilonC = 0;
bHop = 0; % Flag for bound motions on or off
ffo = [0 0.1 0.2 0.3 0.4 0.5]; % filling fraction of obstacles
be = Inf; % binding energy
% Calculated things
numGr = Nr * Nc;
if epsilonC + epsilonR == 0  || epsilonC + epsilonR == 2
  epsilonR = 1;
  epsilonC = 0;
end
% Allocate
energyGrid = zeros( Nr, Nc );
diffMat =  zeros( length(be), length(ffo), numRuns );
% Loop
for ii = 1:length(be)
  beTemp = be(ii);
  for jj = 1:length(ffo)
    nObst = round( numGr * ffo(jj) );
    for kk = 1:numRuns
      [obstGrid, energyGrid] = placeObstacles( nObst, numGr, beTemp );
      diffMat(ii,jj,kk) = genMercSlater( Nr, Nc, numGr, ...
        obstGrid, energyGrid, bHop, epsilonR, epsilonC );
    end
  end
end
% Average
dAve = mean( diffMat, 3 );
dStd = std( diffMat, 0, 3 );

