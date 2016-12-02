function W = construct_weights(regular_grid,extended_grid)

xx = [linspace(1/(regular_grid(1)-extended_grid(1)+1),1,regular_grid(1)-extended_grid(1)),ones(1,regular_grid(2)-regular_grid(1)+1),fliplr(linspace(1/(extended_grid(2)-regular_grid(2)+1),1,extended_grid(2)-regular_grid(2)))];
yy = [linspace(1/(regular_grid(3)-extended_grid(3)+1),1,regular_grid(3)-extended_grid(3)),ones(1,regular_grid(4)-regular_grid(3)+1),fliplr(linspace(1/(extended_grid(4)-regular_grid(4)+1),1,extended_grid(4)-regular_grid(4)))];
zz = [linspace(1/(regular_grid(5)-extended_grid(5)+1),1,regular_grid(5)-extended_grid(5)),ones(1,regular_grid(6)-regular_grid(5)+1),fliplr(linspace(1/(extended_grid(6)-regular_grid(6)+1),1,extended_grid(6)-regular_grid(6)))];

[XX,YY,ZZ] = meshgrid(xx,yy,zz);
W = XX.*YY.*ZZ;