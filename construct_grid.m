function [xx_s,xx_f,yy_s,yy_f,zz_s,zz_f,xx_us,xx_uf,yy_us,yy_uf,zz_us,zz_uf] = construct_grid(grid_size,mot_uf,d1,d2,d3,min_patch_size)

xx_s = 1:grid_size(1):d1;
yy_s = 1:grid_size(2):d2;
zz_s = 1:grid_size(3):d3;

xx_f = [xx_s(2:end)-1,d1];
yy_f = [yy_s(2:end)-1,d2];
zz_f = [zz_s(2:end)-1,d3];

if xx_f(end)-xx_s(end) + 1 < min_patch_size(1) && length(xx_s) > 1; xx_s(end) = []; xx_f(end-1) = []; end
if yy_f(end)-yy_s(end) + 1 < min_patch_size(2) && length(yy_s) > 1; yy_s(end) = []; yy_f(end-1) = []; end
if zz_f(end)-zz_s(end) + 1 < min_patch_size(3) && length(zz_s) > 1; zz_s(end) = []; zz_f(end-1) = []; end

grid_size_us = floor(grid_size./mot_uf);
if mot_uf(1) > 1
    xx_us = 1:grid_size_us(1):d1;
    xx_uf = [xx_us(2:end)-1,d1];
else
    xx_us = xx_s; xx_uf = xx_f;
end
if mot_uf(2) > 1
    yy_us = 1:grid_size_us(2):d2;
    yy_uf = [yy_us(2:end)-1,d2];
else
    yy_us = yy_s; yy_uf = yy_f;
end
if mot_uf(3) > 1
    zz_us = 1:grid_size_us(3):d3;
    zz_uf = [zz_us(2:end)-1,d3];
else
    zz_us = zz_s; zz_uf = zz_f;
end