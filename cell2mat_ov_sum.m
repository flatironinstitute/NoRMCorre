function X = cell2mat_ov_sum(I,xx_s,xx_f,yy_s,yy_f,zz_s,zz_f,overlap,sz,Bs)

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

if nargin < 10 || isempty(Bs)
    Bs = cellfun(@(x) ones(size(x)), I,'un',0);
end
X = zeros([sz,size(I{1,1},length(sz)+1)]);
B = zeros(size(X));
if length(sz) == 2; sz(3) = 1; end

for i = 1:length(xx_f)
    for j = 1:length(yy_f)
        for k = 1:length(zz_f)
            extended_grid = [max(xx_s(i)-overlap(1),1),min(xx_f(i)+overlap(1),sz(1)),max(yy_s(j)-overlap(2),1),min(yy_f(j)+overlap(2),sz(2)),max(zz_s(k)-overlap(3),1),min(zz_f(k)+overlap(3),sz(3))];
            %W = construct_weights([xx_s(i),xx_f(i),yy_s(j),yy_f(j),zz_s(k),zz_f(k)],extended_grid)';            
            Xtemp = zeros(size(I{i,j,k}));
            Btemp = Xtemp;
            ind = ~isnan(I{i,j,k});
            Xtemp(ind) = Bs{i,j,k}(ind).*I{i,j,k}(ind);
            Btemp(ind) = Bs{i,j,k}(ind);
            %X(extended_grid(1):extended_grid(2),extended_grid(3):extended_grid(4),extended_grid(5):extended_grid(6)) = X(extended_grid(1):extended_grid(2),extended_grid(3):extended_grid(4),extended_grid(5):extended_grid(6)) + Bs{i,j,k}.*I{i,j,k}; 
            %B(extended_grid(1):extended_grid(2),extended_grid(3):extended_grid(4),extended_grid(5):extended_grid(6)) = B(extended_grid(1):extended_grid(2),extended_grid(3):extended_grid(4),extended_grid(5):extended_grid(6)) + Bs{i,j,k}.*(I{i,j,k}~=0);
            X(extended_grid(1):extended_grid(2),extended_grid(3):extended_grid(4),extended_grid(5):extended_grid(6)) = X(extended_grid(1):extended_grid(2),extended_grid(3):extended_grid(4),extended_grid(5):extended_grid(6)) + Xtemp; 
            B(extended_grid(1):extended_grid(2),extended_grid(3):extended_grid(4),extended_grid(5):extended_grid(6)) = B(extended_grid(1):extended_grid(2),extended_grid(3):extended_grid(4),extended_grid(5):extended_grid(6)) + Btemp;            
        end
    end
end

X = X./B;
%X(isnan(X))=0;