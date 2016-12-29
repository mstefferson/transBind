%place obstacles on seperate lattice sites
function [obstGrid, energyGrid] = place_obstacles( nObj, numGr, be )
% allocate
energyGrid = zeros( numGr, numGr);
obstGrid = zeros( numGr, numGr);
% indices
inds = randperm( numGr, nObj );
% place obstacles. No overlap
energyGrid( inds ) = be;
obstGrid( inds ) = 1;
end
