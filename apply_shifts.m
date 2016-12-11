function I = apply_shifts(Y,shifts,options)

% apply shifts using fft interpolation

% INPUTS
% Y:                Input data, can be already loaded in memory as a 3D
%                   tensor, a memory mapped file, or a pointer to a tiff stack
% shifts:           calculated shifts
% options:          options structure for motion correction 

% OUTPUTS
% I:                registered data

T = length(shifts);

sizY = size(Y);
if sizY(end) == T
    flag_constant = false;
    nd = length(sizY)-1;   
    sizY = sizY(1:end-1);
else
    flag_constant = true;
    nd = length(sizY);
end

[d1,d2,d3] = size(Y);
if nd == 2; d3 = 1; end

%% precompute some quantities that are used repetitively for template matching and applying shifts

[xx_s,xx_f,yy_s,yy_f,zz_s,zz_f,xx_us,xx_uf,yy_us,yy_uf,zz_us,zz_uf] = construct_grid(options.grid_size,options.mot_uf,d1,d2,d3,options.min_patch_size);
M_fin = mat2cell_ov(zeros(d1,d2,d3),xx_us,xx_uf,yy_us,yy_uf,zz_us,zz_uf,options.overlap_post,[d1,d2,d3]);
temp_cell = M_fin;
Nr = cell(size(temp_cell));
Nc = cell(size(temp_cell));
Np = cell(size(temp_cell));
Bs = cell(size(temp_cell));
for i = 1:length(xx_us)
    for j = 1:length(yy_us)
        for k = 1:length(zz_us)
            [nr,nc,np] = size(temp_cell{i,j,k});
            nr = ifftshift(-fix(nr/2):ceil(nr/2)-1);
            nc = ifftshift(-fix(nc/2):ceil(nc/2)-1);
            np = ifftshift(-fix(np/2):ceil(np/2)-1);
            [Nc{i,j,k},Nr{i,j,k},Np{i,j,k}] = meshgrid(nc,nr,np);
            extended_grid = [max(xx_us(i)-options.overlap_post(1),1),min(xx_uf(i)+options.overlap_post(1),d1),max(yy_us(j)-options.overlap_post(2),1),min(yy_uf(j)+options.overlap_post(2),d2),max(zz_us(k)-options.overlap_post(3),1),min(zz_uf(k)+options.overlap_post(3),d3)];            
            Bs{i,j,k} = permute(construct_weights([xx_us(i),xx_uf(i),yy_us(j),yy_uf(j),zz_us(k),zz_uf(k)],extended_grid),[2,1,3]); 
        end
    end
end
if nd == 2; Np = cellfun(@(x) 0,Nr,'un',0); end

if nd == 2; I = zeros(d1,d2,T); end
if nd == 3; I = zeros(d1,d2,d3,T); end

%%
if flag_constant;  
    Yc = mat2cell_ov(Y,xx_us,xx_uf,yy_us,yy_uf,zz_us,zz_uf,options.overlap_post,[d1,d2,d3]);
    Yfft = cellfun(@(x) fftn(x),Yc,'un',0);
end

for t = 1:T  
    if ~flag_constant
        Yc = mat2cell_ov(Y(:,:,t),xx_us,xx_uf,yy_us,yy_uf,zz_us,zz_uf,options.overlap_post,[d1,d2,d3]);
        Yfft = cellfun(@(x) fftn(x),Yc,'un',0);
    end
    if all(options.mot_uf == 1)
        M_fin = shift_reconstruct2(Yfft{1},shifts(t).shifts,shifts(t).diff,options.us_fac,Nr{1},Nc{1},Np{1},'NaN',0);
        if nd == 2; I(:,:,t) = M_fin; end
        if nd == 3; I(:,:,:,t) = M_fin; end
    else
        shifts_up = shifts(t).shifts_up;
        for i = 1:length(xx_uf)
            for j = 1:length(yy_uf)
                for k = 1:length(zz_uf)                  
                     M_fin{i,j,k} = shift_reconstruct2(Yfft{i,j,k},shifts_up(i,j,k,:),shifts(t).diff(i,j,k),options.us_fac,Nr{i,j,k},Nc{i,j,k},Np{i,j,k},'NaN',0);
                end
            end
        end
          
        gx = max(abs(reshape(diff(shifts_up,[],1),[],1)));
        gy = max(abs(reshape(diff(shifts_up,[],2),[],1)));
        gz = max(abs(reshape(diff(shifts_up,[],3),[],1)));
        flag_interp = max([gx;gy;gz;0])<0.5;      % detect possible smearing

        if flag_interp    
            Mf = cell2mat_ov_sum(M_fin,xx_us,xx_uf,yy_us,yy_uf,zz_us,zz_uf,options.overlap_post,sizY,Bs);
        else            
            Mf = cell2mat_ov(M_fin,xx_us,xx_uf,yy_us,yy_uf,zz_us,zz_uf,options.overlap_post,sizY);
        end    
        if nd == 2; I(:,:,t) = Mf; end
        if nd == 3; I(:,:,:,t) = Mf; end
    end
end
