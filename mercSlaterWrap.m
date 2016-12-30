% Grid
Nc = 4; % number columns
Nr = 2; % number rows
numRuns = 10; % number of runs
% Temporary epsilon
epsilonR = 0; % row bias. if 0, bias along columns
bHop = 0; % Flag for bound motions on or off
ffo = [0]; % filling fraction of obstacles
be = Inf; % binding energy
% Calculated things
numGr = Nr * Nc;
% Allocate
energyGrid = zeros( Nr, Nc );
obstGrid =  zeros( Nr, Nc );
obstGrid(2,1) = 1; obstGrid(2,2) = 1;
energyGrid(2,1) = Inf; energyGrid(2,2) = Inf;
diffMat =  zeros( length(be), length(ffo), numRuns );
diffMatBeta =  zeros( length(be), length(ffo), numRuns );
% Loop
for ii = 1:length(be)
  beTemp = be(ii);
  for jj = 1:length(ffo)
    nObst = round( numGr * ffo(jj) );
    for kk = 1:numRuns
      %[obstGrid, energyGrid] = placeObstacles( nObst, numGr, beTemp );
      diffMat(ii,jj,kk) = genMercSlater( Nr, Nc, numGr, ...
        obstGrid, energyGrid, bHop, epsilonR );
      diffMatBeta(ii,jj,kk) = betaMercSlater( Nr, Nc, numGr, ...
        obstGrid, epsilonR );
    end
  end
end
% Average
dAve = mean( diffMat, 3 );
dStd = std( diffMat, 0, 3 );
% Average
dAveBeta = mean( diffMatBeta, 3 );
dStdBeta = std( diffMatBeta, 0, 3 );

