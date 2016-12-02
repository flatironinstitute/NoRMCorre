function X = cell2mat_ov(I,xx_s,xx_f,yy_s,yy_f,zz_s,zz_f,overlap,sz)

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

X = NaN([sz,size(I{1,1},length(sz)+1)]);
if length(sz) == 2; sz(3) = 1; end

for i = 1:length(xx_f)
    for j = 1:length(yy_f)
        for k = 1:length(zz_f)
            extended_grid = [max(xx_s(i)-overlap(1),1),min(xx_f(i)+overlap(1),sz(1)),max(yy_s(j)-overlap(2),1),min(yy_f(j)+overlap(2),sz(2)),max(zz_s(k)-overlap(3),1),min(zz_f(k)+overlap(3),sz(3))];
            X(xx_s(i):xx_f(i),yy_s(j):yy_f(j),zz_s(k):zz_f(k)) = ... 
                I{i,j,k}(1+(xx_s(i)-extended_grid(1)):end-(extended_grid(2)-xx_f(i)),1+(yy_s(j)-extended_grid(3)):end-(extended_grid(4)-yy_f(j)),1+(zz_s(k)-extended_grid(5)):end-(extended_grid(6)-zz_f(k)));
        end
    end
end