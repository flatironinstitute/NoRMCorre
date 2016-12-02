function [M,shifts,template,shifts_up,M_final] = motion_correction_patches(Y,grid_size,mot_uf,overlap,bin_width,init_batch,us_fac,method)

Y = double(Y);
nY = min(Y(:));
Y = Y - nY;

if nargin < 8 || isempty(method); method = 'mean'; end
if nargin < 7 || isempty(us_fac); us_fac = 10; end
if nargin < 6 || isempty(init_batch); init_batch = 30; end
if nargin < 5 || isempty(bin_width); bin_width = 10; end
if nargin < 4 || isempty(overlap); overlap = 8; end
if nargin < 3 || isempty(mot_uf); mot_uf = 4; end
if nargin < 2 || isempty(grid_size); grid_size = 128; end

[d1,d2,T] = size(Y);

xx = 1:grid_size:d1;
yy = 1:grid_size:d2;

shifts = zeros(length(xx),length(yy),T,2);
shifts_up = zeros(mot_uf*length(xx),mot_uf*length(yy),T,2);
M = zeros(size(Y));
template = zeros(d1,d2);

for i = 1:length(xx)
    for j = 1:length(yy)
        Y_temp = Y(xx(i):min(xx(i)+grid_size-1,d1),yy(j):min(yy(j)+grid_size-1,d2),:);
        [M_temp,shifts_temp,template_temp] = online_motion_correction(Y_temp,bin_width,init_batch,us_fac,method);
        M(xx(i):min(xx(i)+grid_size-1,d1),yy(j):min(yy(j)+grid_size-1,d2),:) = M_temp;
        template(xx(i):min(xx(i)+grid_size-1,d1),yy(j):min(yy(j)+grid_size-1,d2)) = template_temp;
        shifts(i,j,:,:) = shifts_temp;
        disp([i,j])
    end
end

M = M + nY;

for t = 1:T
    for dim = 1:2
        shifts_up(:,:,t,dim) = imresize(shifts(:,:,t,dim),mot_uf);
    end
end

grid_size_fine = floor(grid_size/mot_uf);
xx_fine = 1:grid_size_fine:d1;
yy_fine = 1:grid_size_fine:d2;
M_final = zeros(size(Y));

profile on
for t = 1:T
    Y_temp = Y(:,:,t);
    M_fin = cell(length(xx_fine),length(yy_fine)); %zeros(size(Y_temp));
    for i = 1:length(xx_fine)
        for j = 1:length(yy_fine)
            true_grid = [xx_fine(i),min(xx_fine(i)+grid_size_fine-1,d1),yy_fine(j),min(yy_fine(j)+grid_size_fine-1,d2)];
            extended_grid = [max(xx_fine(i)-overlap,1),min(xx_fine(i)+grid_size_fine+overlap-1,d1),max(yy_fine(j)-overlap,1),min(yy_fine(j)+grid_size_fine+overlap-1,d2)];
            %I_temp = Y(extended_grid(1):extended_grid(2),extended_grid(3):extended_grid(4),t);
            I_temp = Y_temp(extended_grid(1):extended_grid(2),extended_grid(3):extended_grid(4));
            M_temp = shift_reconstruct(I_temp,shifts_up(i,j,t,:),10);
            %M_temp = shift_reconstruct(Y(xx_fine(i):min(xx_fine(i)+grid_size_fine-1,d1),yy_fine(j):min(yy_fine(j)+grid_size_fine-1,d2),t),shifts_up(i,j,t,:),10);
            %M_final(true_grid(1):true_grid(2),true_grid(3):true_grid(4),t) = M_temp(1+(true_grid(1)-extended_grid(1)):end-(extended_grid(2)-true_grid(2)),1+(true_grid(3)-extended_grid(3)):end-(extended_grid(4)-true_grid(4)));
            M_fin{i,j} = M_temp(1+(true_grid(1)-extended_grid(1)):end-(extended_grid(2)-true_grid(2)),1+(true_grid(3)-extended_grid(3)):end-(extended_grid(4)-true_grid(4)));
        end
    end
    M_final(:,:,t) = cell2mat(M_fin);
    if mod(t,100) == 0
        disp(t)
    end
end
profile off;
profile viewer;

M_final = M_final + nY;