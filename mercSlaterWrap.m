% Grid
nC = 3; % number columns
nR = 6; % number rows
numRuns = 1; % number of runs
% Temporary epsilon
epsilonR = 1; % row bias. if 0, bias along columns
bHop = 0; % Flag for bound motions on or off
ffo = [0.2]; % filling fraction of obstacles
be = Inf; % binding energy
% Calculated things
numGr = nR * nC;
% Allocate
myGrid = 1;
if myGrid
  energyGrid = zeros( nR, nC );
  obstGrid =  zeros( nR, nC );
  obstGrid(1,3) = 1; 
  obstGrid(3,:) = 1;
  obstGrid(4,1) = 1;
  obstGrid(5,2:3) = 1;
  energyGrid(obstGrid == 1) = Inf;
end
diffMat =  zeros( length(be), length(ffo), numRuns );
diffMatBeta =  zeros( length(be), length(ffo), numRuns );
% Loop
for ii = 1:length(be)
  beTemp = be(ii);
  for jj = 1:length(ffo)
    nObst = round( numGr * ffo(jj) );
    for kk = 1:numRuns
      if ~myGrid
        [obstGrid, energyGrid] = placeObstacles( nObst, nR, nC, numGr, beTemp );
      end
       diffMat(ii,jj,kk) = genMercSlater( nR, nC, numGr, ...
        obstGrid, energyGrid, bHop, epsilonR );
       diffMatBeta(ii,jj,kk) = betaMercSlater( nR, nC, numGr, ...
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

disp(dAve); disp(dAveBeta);

