function I = mat2cell_ov_2d(X,grid_size,overlap)

% converts a matrix into a cell array with overlapping elements
% INPUTS:
% X:            Input matrix
% grid_size:    size of each element without overlap
% overlap:      amount of overlap

% OUTPUT:
% I:            output cell array

% Written by Eftychios A. Pnevmatikakis, Simons Foundation, 2016

if isscalar(grid_size); grid_size = grid_size*[1,1]; end

[d1,d2,~] = size(X);
xx = 1:grid_size(1):d1;
yy = 1:grid_size(2):d2;
I = cell(length(xx),length(yy));

for i = 1:length(xx)
    for j = 1:length(yy)
        extended_grid = [max(xx(i)-overlap,1),min(xx(i)+grid_size(1)+overlap-1,d1),max(yy(j)-overlap,1),min(yy(j)+grid_size(2)+overlap-1,d2)];
        I{i,j} = X(extended_grid(1):extended_grid(2),extended_grid(3):extended_grid(4),:);
    end
end