function I = mat2cell_ov2(X,xx_s,xx_f,yy_s,yy_f,zz_s,zz_f,sz)

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
for i = 1:length(xx_s)
    for j = 1:length(yy_s)
        for k = 1:length(zz_s)
            if nd == 2
                I{i,j} = X(xx_s(i):xx_f(i),yy_s(j):yy_f(j),:);
            else
                I{i,j,k} = X(xx_s(i):xx_f(i),yy_s(j):yy_f(j),zz_s(k):zz_f(k),:);
            end
        end
    end
end