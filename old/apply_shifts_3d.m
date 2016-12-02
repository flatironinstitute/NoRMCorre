function M_final = apply_shifts_3d(Y,shifts,options)

% online motion correction through DFT subpixel registration
% Based on the dftregistration.m function from Manuel Guizar and Jim Fienup

% INPUTS
% Y:                Input data, can be already loaded in memory as a 3D
%                   tensor, a memory mapped file, or a pointer to a tiff stack
% options:          options structure for motion correction 

% OUTPUTS
% M_final:          motion corrected data
% shifts:           calculated shifts

%% set default parameters if not present

defoptions.memmap = true;                     % save motion corrected file in a mat file
defoptions.grid_size = [64,64,4];             % size of each patch to be corrected separately
defoptions.overlap_post = 8;                  % overlap between subsets within each patch
defoptions.us_fac = 10;                       % upsampling factor for subpixel registration
defoptions.filename = 'motion_corrected.mat'; % filename for motion corrected file

if nargin == 1; options = defoptions; end

if ~isfield(options,'memmap'); options.memmap = defoptions.memmap; end; memmap = options.memmap;
if ~isfield(options,'grid_size'); options.grid_size = defoptions.grid_size; end; grid_size = options.grid_size; 
if ~isfield(options,'overlap_post'); options.overlap_post = defoptions.overlap_post; end; overlap_post = options.overlap_post;
if ~isfield(options,'us_fac'); options.us_fac = defoptions.us_fac; end; us_fac = options.us_fac;
if ~isfield(options,'filename'); options.filename = defoptions.filename; end; filename = options.filename;

if isscalar(grid_size); grid_size = grid_size*[1,1,1]; end

%% determine filetype

if isa(Y,'char')
    [~,~,ext] = fileparts(Y);
    ext = ext(2:end);
    if strcmpi(ext,'tif') || strcmpi(ext,'tiff');
        tiffInfo = imfinfo(Y);
        filetype = 'tif';
        T = length(tiffInfo);
    elseif strcmpi(ext,'mat')
        filetype = 'mem';
        Y = matfile(Y,'Writable',true);
        sizY = size(Y);
        T = sizY(end);
    end    
elseif isobject(Y);
    filetype = 'mem';
    sizY = size(Y,'Y');
    T = sizY(end);
else % array loaded in memory
    filetype = 'mat';
    Y = double(Y);
    sizY = size(Y);
    T = sizY(end);
end

%% read initial batch and compute template

sizY = size(Y,'Y');
[d1,d2,d3,~] = size(Y,'Y');

%% setup grids for patches

xx_fine = 1:grid_size(1):d1;
yy_fine = 1:grid_size(2):d2;
zz_fine = 1:grid_size(3):d3;

if strcmpi(filetype,'mat');
    M_final = zeros(size(Y));
else
    M_final = matfile(filename,'Writable',true);
    if memmap
        M_final.M(d1,d2,d3,T) = uint16(0); %zeros(d1,d2,d3,T);
    end
end
add_value = 0;

for t = 1:T
    switch filetype
        case 'tif'
            Yt = double(imread(Y,'Index',t,'Info',tiffInfo));
        case 'mem'
            Yt = double(Y.Y(:,:,:,t));
        case 'mat'
            Yt = Y(:,:,:,t);
    end        
    Yt = Yt + add_value;

    M_fin = cell(length(xx_fine),length(yy_fine),length(zz_fine)); %zeros(size(Y_temp));
    for i = 1:length(xx_fine)
        for j = 1:length(yy_fine)
            for k = 1:length(zz_fine)
                true_grid = [xx_fine(i),min(xx_fine(i)+grid_size(1)-1,d1),yy_fine(j),min(yy_fine(j)+grid_size(2)-1,d2),zz_fine(k),min(zz_fine(k)+grid_size(3)-1,d3)];
                extended_grid = [max(xx_fine(i)-overlap_post,1),min(xx_fine(i)+grid_size(1)+overlap_post-1,d1),max(yy_fine(j)-overlap_post,1),min(yy_fine(j)+grid_size(2)+overlap_post-1,d2),max(zz_fine(k)-overlap_post,1),min(zz_fine(k)+grid_size(3)+overlap_post-1,d3)];
                I_temp = Yt(extended_grid(1):extended_grid(2),extended_grid(3):extended_grid(4),extended_grid(5):extended_grid(6));
                M_temp = shift_reconstruct_3d(I_temp,shifts(i,j,k,t,:),us_fac,1,1,1);
                M_fin{i,j,k} = M_temp(1+(true_grid(1)-extended_grid(1)):end-(extended_grid(2)-true_grid(2)),1+(true_grid(3)-extended_grid(3)):end-(extended_grid(4)-true_grid(4)),1+(true_grid(5)-extended_grid(5)):end-(extended_grid(6)-true_grid(6)));
            end
        end
    end

    Mf = cell2mat(M_fin) - add_value;
    if strcmpi(filetype,'mat');
        M_final(:,:,:,t) = uint16(Mf);
    elseif memmap
        M_final.M(:,:,:,t) = uint16(Mf);
    end                                
end

if ~strcmpi(filetype,'mat');
    M_final.shifts = shifts;
end

end