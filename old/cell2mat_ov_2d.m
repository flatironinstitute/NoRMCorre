function X = cell2mat_ov_2d(I,grid_size,overlap,d1,d2)

% converts a cell array to a matrix when the cell elements overlap
% INPUTS:
% I:            cell array
% grid_size:    true size of each element
% overlap:      amount of overlap in each direction
% d1:           number of rows of matrix
% d2:           number of columns of matrix

% OUTPUT:
% X:            output matrix

% Written by Eftychios A. Pnevmatikakis, Simons Foundation, 2016

if isscalar(grid_size); grid_size = grid_size*[1,1]; end

xx = 1:grid_size(1):d1;
yy = 1:grid_size(2):d2;
X = zeros(d1,d2,size(I{1,1},3));

for i = 1:length(xx)
    for j = 1:length(yy)
        true_grid = [xx(i),min(xx(i)+grid_size(1)-1,d1),yy(j),min(yy(j)+grid_size(2)-1,d2)];
        extended_grid = [max(xx(i)-overlap,1),min(xx(i)+grid_size(1)+overlap-1,d1),max(yy(j)-overlap,1),min(yy(j)+grid_size(2)+overlap-1,d2)];
        X(true_grid(1):true_grid(2),true_grid(3):true_grid(4)) = I{i,j}(1+(true_grid(1)-extended_grid(1)):end-(extended_grid(2)-true_grid(2)),1+(true_grid(3)-extended_grid(3)):end-(extended_grid(4)-true_grid(4)));
    end
end