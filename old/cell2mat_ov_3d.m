function X = cell2mat_ov_3d(I,grid_size,overlap,d1,d2,d3)

% converts a cell array to a matrix when the cell elements overlap
% INPUTS:
% I:            cell array
% grid_size:    true size of each element
% overlap:      amount of overlap in each direction
% d1:           number of rows of matrix
% d2:           number of columns of matrix
% d3:           number of different planes

% OUTPUT:
% X:            output matrix

% Written by Eftychios A. Pnevmatikakis, Simons Foundation, 2016

if isscalar(grid_size); grid_size = grid_size*[1,1,1]; end

xx = 1:grid_size(1):d1;
yy = 1:grid_size(2):d2;
zz = 1:grid_size(3):d3;
X = zeros(d1,d2,d3,size(I{1,1},4));

for i = 1:length(xx)
    for j = 1:length(yy)
        for k = 1:length(zz)
            true_grid = [xx(i),min(xx(i)+grid_size(1)-1,d1),yy(j),min(yy(j)+grid_size(2)-1,d2),zz(k),min(zz(k)+grid_size(3)-1,d3)];
            extended_grid = [max(xx(i)-overlap,1),min(xx(i)+grid_size(1)+overlap-1,d1),max(yy(j)-overlap,1),min(yy(j)+grid_size(2)+overlap-1,d2),max(zz(k)-overlap,1),min(zz(k)+grid_size(3)+overlap-1,d3)];
            X(true_grid(1):true_grid(2),true_grid(3):true_grid(4),true_grid(5):true_grid(6)) ...
                = I{i,j,k}(1+(true_grid(1)-extended_grid(1)):end-(extended_grid(2)-true_grid(2)),1+(true_grid(3)-extended_grid(3)):end-(extended_grid(4)-true_grid(4)),1+(true_grid(5)-extended_grid(5)):end-(extended_grid(6)-true_grid(6)));
        end
    end
end