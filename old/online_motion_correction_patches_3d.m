function [M_final,shifts,template,shifts_up] = online_motion_correction_patches_3d(Y,options,template_in)

% online motion correction through DFT subpixel registration
% Based on the dftregistration.m function from Manuel Guizar and Jim Fienup

% INPUTS
% Y:                Input data, can be already loaded in memory as a 3D
%                   tensor, a memory mapped file, or a pointer to a tiff stack
% options:          options structure for motion correction 

% OUTPUTS
% M_final:          motion corrected data
% shifts_up:        upsampled shifts
% shifts:           originally calculated shifts
% template:         calculated template


%% set default parameters if not present

defoptions.memmap = true;                     % save motion corrected file in a mat file
defoptions.grid_size = [64,64,4];             % size of each patch to be corrected separately
defoptions.mot_uf = [1,1,1];                  % upsampling factor within each patch
defoptions.overlap_pre = 8;                   % overlap between subsets within each patch
defoptions.overlap_post = 8;                  % overlap between subsets within each patch
defoptions.bin_width = 10;                    % width of buffer for computing the moving template
defoptions.buffer_width = 100;                % number of local means to keep in memory
defoptions.init_batch = 30;                   % length of initial batch
defoptions.us_fac = 10;                       % upsampling factor for subpixel registration
defoptions.method = {'median';'mean'};        % method for averaging the template
defoptions.plot_flag = false;                 % flag for plotting results in real time
defoptions.filename = 'motion_corrected.mat'; % filename for motion corrected file
defoptions.use_parallel = false;              % use parfor when breaking each frame into patches
defoptions.iter = 1;                          % number of passes over the data
%defoptions.write_tiff = false;                % save output as a tiff stack
%defoptions.out_name = 'motion_corrected.tif'; % name for output file name

if nargin == 1; options = defoptions; end

if ~isfield(options,'memmap'); options.memmap = defoptions.memmap; end; memmap = options.memmap;
if ~isfield(options,'grid_size'); options.grid_size = defoptions.grid_size; end; grid_size = options.grid_size; 
if ~isfield(options,'mot_uf'); options.mot_uf = defoptions.mot_uf; end; mot_uf = options.mot_uf;
if ~isfield(options,'overlap_pre'); options.overlap_pre = defoptions.overlap_pre; end; overlap_pre = options.overlap_pre;
if ~isfield(options,'overlap_post'); options.overlap_post = defoptions.overlap_post; end; overlap_post = options.overlap_post;
if ~isfield(options,'bin_width'); options.bin_width = defoptions.bin_width; end; bin_width = options.bin_width;
if ~isfield(options,'buffer_width'); options.buffer_width = defoptions.buffer_width; end; buffer_width = options.buffer_width;
if ~isfield(options,'init_batch'); options.init_batch = defoptions.init_batch; end; init_batch = options.init_batch;
if ~isfield(options,'us_fac'); options.us_fac = defoptions.us_fac; end; us_fac = options.us_fac;
if ~isfield(options,'method'); options.method = defoptions.method; end; method = options.method;
if ~isfield(options,'filename'); options.filename = defoptions.filename; end; filename = options.filename;
if ~isfield(options,'use_parallel'); options.use_parallel = defoptions.use_parallel; end; use_parallel = options.use_parallel;
if ~isfield(options,'iter'); options.iter = defoptions.iter; end; iter = options.iter;
%if ~isfield(options,'write_tiff'); options.write_tiff = defoptions.write_tiff; end; write_tiff = options.write_tiff;
%if ~isfield(options,'out_name'); options.out_name = defoptions.out_name; end; out_name = options.out_name;

if isscalar(grid_size); grid_size = grid_size*[1,1,1]; end
if isscalar(mot_uf); mot_uf = mot_uf*[1,1,1]; end
if ~isfield(options,'max_shift'); options.max_shift = grid_size./mot_uf; end; max_shift = options.max_shift;

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

switch filetype
    case 'tif'
        Y_temp = double(bigread2(Y,1,init_batch));
    case 'mem'
        Y_temp = double(Y.Y(:,:,:,1:init_batch));
    case 'mat'
        Y_temp = Y(:,:,:,1:init_batch);
end

add_value = 0; %mean(Y_temp(:));

if nargin < 3 || isempty(template_in);
    template_in = median(Y_temp,4)+add_value;
end

[d1,d2,d3,~] = size(Y_temp);

%% setup grids for patches

xx = 1:grid_size(1):d1;
yy = 1:grid_size(2):d2;
zz = 1:grid_size(3):d3;

shifts = zeros(length(xx),length(yy),length(zz),T,3);
shifts_up = zeros(mot_uf(1)*length(xx),mot_uf(2)*length(yy),mot_uf(3)*length(zz),T,3);

grid_size_fine = floor(grid_size./mot_uf);
xx_fine = 1:grid_size_fine(1):d1;
yy_fine = 1:grid_size_fine(2):d2;
zz_fine = 1:grid_size_fine(3):d3;

%% precompute some quantities that are used repetitively for template matching
temp_cell = mat2cell_ov(template_in,grid_size_fine,overlap_post);
Nr = cell(size(temp_cell));
Nc = cell(size(temp_cell));
Np = cell(size(temp_cell));
for i = 1:length(xx_fine)
    for j = 1:length(yy_fine)
        for k = 1:length(zz_fine)
            [nr,nc,np] = size(temp_cell{i,j});
            nr = ifftshift(-fix(nr/2):ceil(nr/2)-1);
            nc = ifftshift(-fix(nc/2):ceil(nc/2)-1);
            np = ifftshift(-fix(np/2):ceil(np/2)-1);
            [Nc{i,j},Nr{i,j},Np{i,j}] = meshgrid(nc,nr,np);
        end
    end
end

%%

template = mat2cell_ov_3d(template_in,grid_size,overlap_pre);
fftTemp = cellfun(@fftn,template,'un',0);
buffer = mat2cell_ov_3d(zeros(d1,d2,d3,bin_width),grid_size,overlap_pre);
if strcmpi(method{2},'median')
    buffer_med = zeros(d1,d2,d3,buffer_width);
end

if strcmpi(filetype,'mat');
    M_final = zeros(size(Y));
else
    M_final = matfile(filename,'Writable',true);
    if memmap
        M_final.M(d1,d2,d3,T) = uint16(0);
    end
end

cnt_buf = 0;

for it = 1:iter
    for t = 1:T
        if mod(t,bin_width) == 0
            disp(t)
            cnt_buf = cnt_buf + 1;                
            if strcmpi(method{2},'mean')
                new_temp = cellfun(@(x) mean(x,4), buffer, 'UniformOutput',false);
            elseif strcmpi(method{2},'median');
                new_temp = cellfun(@(x) median(x,4), buffer, 'UniformOutput', false);
            end
            if strcmpi(method{1},'mean')
                cnt = t/bin_width + 1;
                template = cellfun(@plus, cellfun(@(x) x*(cnt-1)/cnt, template,'un',0), cellfun(@(x) x*1/cnt, new_temp,'un',0), 'un',0);
            elseif strcmpi(method{1},'median');
                if cnt_buf <= buffer_width
                    buffer_med(:,:,:,cnt_buf) = cell2mat_ov_3d(new_temp,grid_size,overlap_pre,d1,d2,d3);
                else
                    buffer_med = circshift(buffer_med,[0,0,0,-1]);
                    buffer_med(:,:,:,buffer_width) = cell2mat_ov_3d(new_temp,grid_size,overlap_pre,d1,d2,d3);
                end
                template = mat2cell_ov_3d(median(buffer_med,4),grid_size,overlap_pre);
            end
            fftTemp = cellfun(@fftn, template, 'un',0);
        end

        switch filetype
            case 'tif'
                Yt = double(imread(Y,'Index',t,'Info',tiffInfo));
            case 'mem'
                Yt = double(Y.Y(:,:,:,t));
            case 'mat'
                Yt = Y(:,:,:,t);
        end        
        Yt = Yt + add_value;
        ind = rem(t,bin_width) + bin_width*(rem(t,bin_width)==0);
        Yc = mat2cell_ov_3d(Yt,grid_size,overlap_pre);
        fftY = cellfun(@fftn, Yc, 'un',0);
        if ~use_parallel
            for i = 1:length(xx)
                for j = 1:length(yy)                
                    for k = 1:length(zz)
                        [output,Greg] = dftregistration_max_3d(fftTemp{i,j,k},fftY{i,j,k},us_fac,max_shift);
                        M_temp = abs(ifftn(Greg));
                        buffer{i,j,k}(:,:,:,ind) = M_temp;
                        shifts(i,j,k,t,1) = output(3);
                        shifts(i,j,k,t,2) = output(4);
                        shifts(i,j,k,t,3) = output(5);
                    end
                end
            end
        else
            Mt2 = cell(length(xx)*length(yy)*length(zz),1);
            shifts_temp = cell(length(xx)*length(yy)*length(zz),1);
            parfor ii = 1:length(xx)*length(yy)*length(zz)
                [i,j,k] = ind2sub([length(xx),length(yy),length(zz)],ii)
                [output,Greg] = dftregistration_max_3d(fftTemp{i,j,k},fftY{i,j,k},us_fac,max_shift);
                M_temp = abs(ifftn(Greg));
                Mt2{ii} = M_temp;
                shifts_temp{ii} = output(3:5);
            end
            for ii = 1:length(xx)*length(yy)*length(zz)
                 [i,j,k] = ind2sub([length(xx),length(yy),length(zz)],ii);
                 buffer{i,j,k}(:,:,ind) = Mt2{ii};
                 shifts(i,j,k,t,1) = shifts_temp{ii}(1);
                 shifts(i,j,k,t,2) = shifts_temp{ii}(2);
                 shifts(i,j,k,t,3) = shifts_temp{ii}(3);
            end
        end       
        shifts_p = imresize(shifts(:,:,:,t,1),[mot_uf(1)*length(xx),mot_uf(2)*length(yy)]);
        if mot_uf(3) > 1
            shifts_p = reshape(shifts_p,[],size(shifts_p,3));
            shifts_p = imresize(shifts_p,size(shifts_p).*[1,mot_uf(3)]);
            shifts_p = reshape(shifts_p,mot_uf(1)*length(xx),mot_uf(2)*length(yy),mot_uf(3)*length(zz));
        end
        shifts_up(:,:,:,t,1) = shifts_p;
        shifts_p = imresize(shifts(:,:,:,t,2),[mot_uf(1)*length(xx),mot_uf(2)*length(yy)]);
        if mot_uf(3) > 1
            shifts_p = reshape(shifts_p,[],size(shifts_p,3));
            shifts_p = imresize(shifts_p,size(shifts_p).*[1,mot_uf(3)]);
            shifts_p = reshape(shifts_p,mot_uf(1)*length(xx),mot_uf(2)*length(yy),mot_uf(3)*length(zz));
            shifts_up(:,:,:,t,2) = shifts_p;
        end
%       shifts_up(:,:,t,2) = imresize(shifts(:,:,t,2),mot_uf);

        M_fin = cell(length(xx_fine),length(yy_fine),length(zz_fine)); %zeros(size(Y_temp));
        for i = 1:length(xx_fine)
            for j = 1:length(yy_fine)
                for k = 1:length(zz_fine)
                    true_grid = [xx_fine(i),min(xx_fine(i)+grid_size_fine(1)-1,d1),yy_fine(j),min(yy_fine(j)+grid_size_fine(2)-1,d2),zz_fine(k),min(zz_fine(k)+grid_size_fine(3)-1,d3)];
                    extended_grid = [max(xx_fine(i)-overlap_post,1),min(xx_fine(i)+grid_size_fine(1)+overlap_post-1,d1),max(yy_fine(j)-overlap_post,1),min(yy_fine(j)+grid_size_fine(2)+overlap_post-1,d2),max(zz_fine(k)-overlap_post,1),min(zz_fine(k)+grid_size_fine(3)+overlap_post-1,d3)];
                    I_temp = Yt(extended_grid(1):extended_grid(2),extended_grid(3):extended_grid(4),extended_grid(5):extended_grid(6));
                    M_temp = shift_reconstruct_3d(I_temp,shifts_up(i,j,k,t,:),us_fac,1,1,1,Nr{i,j},Nc{i,j},Np{i,j});
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
        disp(t)
    end
end
template = cellfun(@(x) x - add_value,template,'un',0);
template = cell2mat_ov_3d(template,grid_size,overlap_pre,d1,d2,d3);
if ~strcmpi(filetype,'mat');
    M_final.shifts = shifts;
    M_final.shifts_up = shifts_up; 
    M_final.template = template;
end

end