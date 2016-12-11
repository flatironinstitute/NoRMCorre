function M = apply_shifts_dev(Y,shifts,overlap,xx_us,xx_uf,yy_us,yy_uf,zz_us,zz_uf,us_fac,Nr,Nc,Np,Bs)

% apply shifts using an fft method

% INPUTS
% Y:                Input data, can be already loaded in memory as a 3D
%                   tensor, a memory mapped file, or a pointer to a tiff stack
% options:          options structure for motion correction 

% OUTPUTS
% M_final:          motion corrected data
% shifts:           calculated shifts

%% read initial batch and compute template

if nargin < 10 || isempty(us_fac); us_fac = 2; end
if nargin < 9 || isempty(zz_uf); zz_uf = 1; end
if nargin < 8 || isempty(zz_us); zz_us = 1; end

[d1,d2,d3] = size(Y);
sizY = [d1,d2,d3];

M_fin = cell(length(xx_uf),length(yy_uf),length(zz_uf));

for i = 1:length(xx_uf)
    for j = 1:length(yy_uf)
        for k = 1:length(zz_uf)
            extended_grid = [max(xx_us(i)-overlap(1),1),min(xx_uf(i)+overlap(1),d1),max(yy_us(j)-overlap(2),1),min(yy_uf(j)+overlap(2),d2),max(zz_us(k)-overlap(3),1),min(zz_uf(k)+overlap(3),d3)];
            I_temp = Y(extended_grid(1):extended_grid(2),extended_grid(3):extended_grid(4),extended_grid(5):extended_grid(6));
            M_fin{i,j,k} = shift_reconstruct(I_temp,shifts(i,j,k,:),us_fac,Nr{i,j,k},Nc{i,j,k},Np{i,j,k});
        end
    end
end

gx = max(abs(reshape(diff(shifts,[],1),[],1)));
gy = max(abs(reshape(diff(shifts,[],2),[],1)));
gz = max(abs(reshape(diff(shifts,[],3),[],1)));
flag_interp = max([gx;gy;gz;0])<0.5;      % detect possible smearing

if flag_interp    
    M = cell2mat_ov_sum(M_fin,xx_us,xx_uf,yy_us,yy_uf,zz_us,zz_uf,overlap,sizY,Bs);
else            
    M = cell2mat_ov(M_fin,xx_us,xx_uf,yy_us,yy_uf,zz_us,zz_uf,overlap,sizY);
end