%place obstacles on seperate lattice sites
function [obstGrid, energyGrid] = placeObstacles( nObj, nR, nC, numGr, be )
% allocate
energyGrid = zeros( nR, nC);
obstGrid = zeros( nR, nC);
% indices
inds = randperm( numGr, nObj );
% place obstacles. No overlap
energyGrid( inds ) = be;
obstGrid( inds ) = 1;
end
