function X = cell2mat_ov2(I,xx_s,xx_f,yy_s,yy_f,zz_s,zz_f,overlap,sz)

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

X = zeros([sz,size(I{1,1},length(sz)+1)]);

for i = 1:length(xx_f)
    for j = 1:length(yy_f)
        for k = 1:length(zz_f)
            regular_grid = [xx_s(i)+overlap(1)*(i>1),xx_f(i)-overlap(1)*(i<length(xx_f)),yy_s(i)+overlap(2)*(j>1),jj_f(j)-overlap(2)*(j<length(yy_f)),zz_s(k)+overlap(3)*(k>1),zz_f(k)-overlap(3)*(k<length(zz_f))];
            %extended_grid = [max(xx_s(i)-overlap(1),1),min(xx_f(i)+overlap(1),sz(1)),max(yy_s(j)-overlap(2),1),min(yy_f(j)+overlap(2),sz(2)),max(zz_s(k)-overlap(3),1),min(zz_f(k)+overlap(3),sz(3))];
            X(regular_grid(1):regular_grid(2),regular_grid(3):regular_grid(4),regular_grid(5):regular_grid(6)) = ... 
                I{i,j,k}(1+(regular_grid(1)-xx_s(i)):end-(xx_f(i)-regular_grid(2)),1+(regular_grid(3)-yy_s(j)):end-(yy_f(j)-regular_grid(4)),1+(regular_grid(5)-zz_s(k)):end-(zz_f(k)-regular_grid(6)));
        end
    end
end