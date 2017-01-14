function I = apply_shifts(Y,shifts,options,td1,td2,td3)

% apply shifts using fft interpolation

% INPUTS
% Y:                Input data, can be already loaded in memory as a 3D
%                   tensor, a memory mapped file, or a pointer to a tiff stack
% shifts:           calculated shifts
% options:          options structure for motion correction 
% td1,td2,td3:      extend patches on the boundaries by that much        

% OUTPUTS
% I:                registered data

if nargin < 6 || isempty(td3); td3 = 0; end
if nargin < 5 || isempty(td2); td2 = 0; end
if nargin < 4 || isempty(td1); td1 = 0; end

T = length(shifts);

sizY = size(Y);
if sizY(end) == T && T > 1
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

[xx_s,xx_f,yy_s,yy_f,zz_s,zz_f,xx_us,xx_uf,yy_us,yy_uf,zz_us,zz_uf] = construct_grid(options.grid_size,options.mot_uf,options.d1,options.d2,options.d3,options.min_patch_size);
xx_us = xx_us + td1; xx_us(1) = 1;
yy_us = yy_us + td2; yy_us(1) = 1;
zz_us = zz_us + td3; zz_us(1) = 1;
xx_uf = xx_uf + td1; xx_uf(end) = d1;
yy_uf = yy_uf + td2; yy_uf(end) = d2;
zz_uf = zz_uf + td3; zz_uf(end) = d3;

M_fin = mat2cell_ov(zeros(d1,d2,d3,'single'),xx_us,xx_uf,yy_us,yy_uf,zz_us,zz_uf,options.overlap_post,[d1,d2,d3]);
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
bin_width = options.bin_width;
for t = 1:bin_width:T
    if nd == 2; Ytm = single(Y(:,:,t:min(t+bin_width-1,T))); end
    if nd == 3; Ytm = single(Y(:,:,:,t:min(t+bin_width-1,T))); end
    if nd == 2; Ytc = mat2cell(Ytm,d1,d2,ones(1,size(Ytm,3))); end
    if nd == 3; Ytc = mat2cell(Ytm,d1,d2,d3,ones(1,size(Ytm,4))); end
    
    Mf = cell(size(Ytc));
    lY = length(Ytc);
    shifts_temp = shifts(t:t+lY-1);
    
    %buffer = cell(length(xx_us),length(yy_us),length(zz_us),size(Ytm,ndims(Ytm)));
    %shifts = struct('shifts',cell(lY,1),'shifts_up',cell(lY,1),'diff',cell(lY,1));
    for ii = 1:lY        
        shifts_temp(ii).diff(:) = 0;
        if ~flag_constant            
            Yc = mat2cell_ov(Ytc{ii},xx_us,xx_uf,yy_us,yy_uf,zz_us,zz_uf,options.overlap_post,[d1,d2,d3]);
            Yfft = cellfun(@(x) fftn(x),Yc,'un',0);
        end
        if all(options.mot_uf == 1)
            M_fin = shift_reconstruct2(Yfft{1},shifts_temp(ii).shifts,shifts_temp(ii).diff,options.us_fac,Nr{1},Nc{1},Np{1},'NaN',0);
            if nd == 2; Mf{ii} = M_fin; end
            if nd == 3; Mf{ii} = M_fin; end
        else
            shifts_up = shifts_temp(ii).shifts_up;
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
                Mf{ii} = cell2mat_ov_sum(M_fin,xx_us,xx_uf,yy_us,yy_uf,zz_us,zz_uf,options.overlap_post,sizY,Bs);
            else            
                Mf{ii} = cell2mat_ov(M_fin,xx_us,xx_uf,yy_us,yy_uf,zz_us,zz_uf,options.overlap_post,sizY);
            end    
        end
    end
    Mf = cell2mat(Mf);
    if nd == 2; I(:,:,t:min(t+bin_width-1,T)) = Mf; end
    if nd == 3; I(:,:,:,t:min(t+bin_width-1,T)) = Mf; end    
end
