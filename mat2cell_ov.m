function I = mat2cell_ov(X,xx_s,xx_f,yy_s,yy_f,zz_s,zz_f,overlap,sz)

% converts a matrix into a cell array with overlapping elements
% INPUTS:
% X:            Input matrix
% grid_size:    size of each element without overlap
% overlap:      amount of overlap
% sz:           spatial size of X

% OUTPUT:
% I:            output cell array

% Written by Eftychios A. Pnevmatikakis, Simons Foundation, 2016

I = cell(length(xx_s),length(yy_s),length(zz_s));
nd = length(sz);
if nd == 2; sz(3) = 1; end
for i = 1:length(xx_s)
    for j = 1:length(yy_s)
        for k = 1:length(zz_s)
            extended_grid = [max(xx_s(i)-overlap(1),1),min(xx_f(i)+overlap(1),sz(1)),max(yy_s(j)-overlap(2),1),min(yy_f(j)+overlap(2),sz(2)),max(zz_s(k)-overlap(3),1),min(zz_f(k)+overlap(3),sz(3))];            
            if nd == 2
                I{i,j} = X(extended_grid(1):extended_grid(2),extended_grid(3):extended_grid(4),:);
            else
                I{i,j,k} = X(extended_grid(1):extended_grid(2),extended_grid(3):extended_grid(4),extended_grid(5):extended_grid(6),:);
            end
        end
    end
end