function [M_final,shifts,template,shifts_up] = online_motion_correction_patches_2d(Y,options,template)

% online motion correction through DFT subpixel registration
% Based on the dftregistration.m function from Manuel Guizar and Jim Fienup

% INPUTS
% Y:                Input data, can be already loaded in memory as a 3D
%                   tensor, a memory mapped file, or a pointer to a tiff stack
% options:          options structure for motion correction 
% template:         provide template (optional)

% OUTPUTS
% M_final:          motion corrected data
% shifts_up:        upsampled shifts
% shifts:           originally calculated shifts
% template:         calculated template

%% set default parameters if not present

defoptions.memmap = true;                     % save motion corrected file in a mat file
defoptions.grid_size = [128,128];             % size of each patch to be corrected separately
defoptions.mot_uf = 4;                        % upsampling factor within each patch
defoptions.overlap_pre = 8;                   % overlap between subsets within each patch
defoptions.overlap_post = 8;                  % overlap between subsets within each patch
defoptions.bin_width = 10;                    % width of buffer for computing the moving template
defoptions.buffer_width = 100;                % number of local means to keep in memory
defoptions.init_batch = 30;                   % length of initial batch
defoptions.us_fac = 5;                        % upsampling factor for subpixel registration
defoptions.method = {'median';'mean'};        % method for averaging the template
defoptions.plot_flag = false;                 % flag for plotting results in real time
defoptions.filename = 'motion_corrected.mat'; % filename for motion corrected mat file
defoptions.use_parallel = false;              % use parfor when breaking each frame into patches
defoptions.make_avi = false;                  % flag for making movie
defoptions.name = 'motion_corrected.avi';     % name of saved movie
defoptions.fr = 15;                           % frame rate for saved movie
defoptions.iter = 1;                          % number of passes over the data
%defoptions.write_tiff = false;                % save output as a tiff stack
%defoptions.out_name = 'motion_corrected.tif'; % name for output file name

if nargin == 1 || isempty(options); options = defoptions; end

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
if ~isfield(options,'plot_flag'); options.plot_flag = defoptions.plot_flag; end; plot_flag = options.plot_flag;
if ~isfield(options,'filename'); options.filename = defoptions.filename; end; filename = options.filename;
if ~isfield(options,'use_parallel'); options.use_parallel = defoptions.use_parallel; end; use_parallel = options.use_parallel;
if ~isfield(options,'make_avi'); options.make_avi = defoptions.make_avi; end; make_avi = options.make_avi;
if ~isfield(options,'name'); options.name = defoptions.name; end; name = options.name;
if ~isfield(options,'fr'); options.fr = defoptions.fr; end; fr = options.fr;
if ~isfield(options,'iter'); options.iter = defoptions.iter; end; iter = options.iter;
%if ~isfield(options,'write_tiff'); options.write_tiff = defoptions.write_tiff; end; write_tiff = options.write_tiff;
%if ~isfield(options,'out_name'); options.out_name = defoptions.out_name; end; out_name = options.out_name;
if ~isfield(options,'max_shift'); options.max_shift = grid_size/mot_uf; end; max_shift = options.max_shift;

if isscalar(grid_size); grid_size = grid_size*[1,1]; end

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
    elseif strcmpi(ext,'hdf5') || strcmpi(ext,'h5');
        filetype = 'hdf5';
        fileinfo = hdf5info(Y);
        sizY = fileinfo.GroupHierarchy.Datasets.Dims;
        T = sizY(end);
    end    
elseif isobject(Y);
    filetype = 'mem';
    sizY = size(Y);
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
    case 'hdf5'
        Y_temp = double(bigread2(Y,1,init_batch));        
    case 'mem'
        Y_temp = double(Y.Y(:,:,1:init_batch));
    case 'mat'
        Y_temp = Y(:,:,1:init_batch);
end

add_value = mean(Y_temp(:));
if nargin < 3 || isempty(template)
    template_in = median(Y_temp,3)+add_value;
else
    template_in = template + add_value;
end

[d1,d2,~] = size(Y_temp);

%% setup grids for patches

xx = 1:grid_size(1):d1;
yy = 1:grid_size(2):d2;

if length(xx) == 1 && length(yy) == 1; mot_uf = 1; end

shifts = zeros(length(xx),length(yy),T,2);
shifts_up = zeros(mot_uf*length(xx),mot_uf*length(yy),T,2);

grid_size_fine = floor(grid_size/mot_uf);
xx_fine = 1:grid_size_fine(1):d1;
yy_fine = 1:grid_size_fine(2):d2;

%% precompute some quantities that are used repetitively for template matching
temp_cell = mat2cell_ov(template_in,grid_size_fine,overlap_post);
Nr = cell(size(temp_cell));
Nc = cell(size(temp_cell));
for i = 1:length(xx_fine)
    for j = 1:length(yy_fine)
        [nr,nc] = size(temp_cell{i,j});
        nr = ifftshift(-fix(nr/2):ceil(nr/2)-1);
        nc = ifftshift(-fix(nc/2):ceil(nc/2)-1);
        [Nc{i,j},Nr{i,j}] = meshgrid(nc,nr);
    end
end
Np = cellfun(@(x) 0,Nr,'un',0);
%%
maxNumCompThreads(2);
template = mat2cell_ov(template_in,grid_size,overlap_pre);
fftTemp = cellfun(@fft2,template,'un',0);
buffer_pre = 1;
if buffer_pre
    buffer = mat2cell_ov(zeros(d1,d2,bin_width),grid_size,overlap_pre);
else
    buffer = mat2cell_ov(zeros(d1,d2,bin_width),grid_size_fine,overlap_post);
end
if strcmpi(method{2},'median')
    buffer_med = zeros(d1,d2,buffer_width);
end

if strcmpi(filetype,'mat');
    M_final = zeros(size(Y));
else
    M_final = matfile(filename,'Writable',true);
    if memmap
        M_final.M = zeros(d1,d2,T);
    end
end

if plot_flag
    if make_avi
        vidObj = VideoWriter(name);
        set(vidObj,'FrameRate',fr);
        open(vidObj);
    end
    if strcmpi(filetype,'mat')
        nnY = quantile(Y(:),0.005);
        mmY = quantile(Y(:),0.995);
    else
        nnY = quantile(Y_temp(:),0.005);
        mmY = quantile(Y_temp(:),0.995);
    end
    fig = figure;
        screensize = get(0,'Screensize' );
        fac = min(min((screensize(3:4)-200)./[d2,d1]),10);
        set(gcf, 'PaperUnits', 'points', 'Units', 'points');
        set(gcf, 'Position', round([100 100 fac*d2 fac*d1]));
end
cnt_buf = 0;
fprintf('Template initialization complete. \n')
for it = 1:iter
    if it < iter; plot_flag = 0; else plot_flag = options.plot_flag; end
    for t = 1:T
        if mod(t,bin_width) == 0
            disp(t)
            cnt_buf = cnt_buf + 1;                
            if strcmpi(method{2},'mean')
                new_temp = cellfun(@(x) mean(x,3), buffer, 'UniformOutput',false);
            elseif strcmpi(method{2},'median');
                new_temp = cellfun(@(x) median(x,3), buffer, 'UniformOutput', false);
            end
            if strcmpi(method{1},'mean')
                cnt = t/bin_width + 1;
                template = cellfun(@plus, cellfun(@(x) x*(cnt-1)/cnt, template,'un',0), cellfun(@(x) x*1/cnt, new_temp,'un',0), 'un',0);
            elseif strcmpi(method{1},'median');
                if cnt_buf <= buffer_width
                    buffer_med(:,:,cnt_buf) = cell2mat_ov(new_temp,grid_size*buffer_pre + grid_size_fine*(1-buffer_pre),overlap_pre*buffer_pre + overlap_post*(1-buffer_pre),d1,d2);
                else
                    buffer_med = circshift(buffer_med,[0,0,-1]);
                    buffer_med(:,:,buffer_width) = cell2mat_ov(new_temp,grid_size*buffer_pre + grid_size_fine*(1-buffer_pre),overlap_pre*buffer_pre + overlap_post*(1-buffer_pre),d1,d2);
                end
                template = mat2cell_ov(median(buffer_med,3),grid_size,overlap_pre);
            end
            fftTemp = cellfun(@fft2, template, 'un',0);
        end

        switch filetype
            case 'tif'
                Yt = double(imread(Y,'Index',t,'Info',tiffInfo));
            case 'hdf5'
                Yt = double(h5read(Y,'/mov',[ones(1,length(sizY)-1),t],[sizY(1:end-1),1]));
            case 'mem'
                Yt = double(Y.Y(:,:,t));
            case 'mat'
                Yt = Y(:,:,t);
        end        
        Yt = Yt + add_value;
        ind = rem(t,bin_width) + bin_width*(rem(t,bin_width)==0);
        Yc = mat2cell_ov(Yt,grid_size,overlap_pre);
        fftY = cellfun(@fft2, Yc, 'un',0);
        
        M_fin = cell(length(xx_fine),length(yy_fine)); %zeros(size(Y_temp));
        if ~use_parallel
            for i = 1:length(xx)
                for j = 1:length(yy)                                                        
                    if buffer_pre
                        [output,Greg] = dftregistration_max(fftTemp{i,j},fftY{i,j},us_fac,max_shift);
                        M_temp = abs(ifft2(Greg));
                        buffer{i,j}(:,:,ind) = M_temp;
                        if mot_uf == 1
                            M_fin{i,j} = M_temp;
                        end
                    else
                        output = dftregistration_max(fftTemp{i,j},fftY{i,j},us_fac,max_shift);
                    end
                    shifts(i,j,t,1) = output(3);
                    shifts(i,j,t,2) = output(4);
                end
            end
        else
            Mt2 = cell(length(xx)*length(yy),1);
            shifts_temp = cell(length(xx)*length(yy),1);
            parfor ii = 1:length(xx)*length(yy)
                [i,j] = ind2sub([length(xx),length(yy)],ii)
                               
                if buffer_pre
                    [output,Greg] = dftregistration_max(fftTemp{i,j},fftY{i,j},us_fac,max_shift); 
                    M_temp = abs(ifft2(Greg));
                    Mt2{ii} = M_temp;
                else
                    output = dftregistration_max(fftTemp{i,j},fftY{i,j},us_fac,max_shift);
                end
                shifts_temp{ii} = output(3:4);
            end
            for ii = 1:length(xx)*length(yy)
                 [i,j] = ind2sub([length(xx),length(yy)],ii);
                 if buffer_pre
                    buffer{i,j}(:,:,ind) = Mt2{ii};
                 end
                 if mot_uf == 1
                     M_fin{i,j} = Mt2{ii};
                 end
                 shifts(i,j,t,1) = shifts_temp{ii}(1);
                 shifts(i,j,t,2) = shifts_temp{ii}(2);
            end
        end       
        
        if mot_uf > 1
            shifts_up(:,:,t,1) = imresize(shifts(:,:,t,1),mot_uf);
            shifts_up(:,:,t,2) = imresize(shifts(:,:,t,2),mot_uf);

            for i = 1:length(xx_fine)
                for j = 1:length(yy_fine)
                    %true_grid = [xx_fine(i),min(xx_fine(i)+grid_size_fine(1)-1,d1),yy_fine(j),min(yy_fine(j)+grid_size_fine(2)-1,d2)];
                    extended_grid = [max(xx_fine(i)-overlap_post,1),min(xx_fine(i)+grid_size_fine(1)+overlap_post-1,d1),max(yy_fine(j)-overlap_post,1),min(yy_fine(j)+grid_size_fine(2)+overlap_post-1,d2)];
                    I_temp = Yt(extended_grid(1):extended_grid(2),extended_grid(3):extended_grid(4));
                    M_temp = shift_reconstruct_3d(I_temp,shifts_up(i,j,t,:),us_fac,1,1,1,Nr{i,j},Nc{i,j},Np{i,j});
                    M_fin{i,j} = M_temp; %(1+(true_grid(1)-extended_grid(1)):end-(extended_grid(2)-true_grid(2)),1+(true_grid(3)-extended_grid(3)):end-(extended_grid(4)-true_grid(4)));
                    if ~buffer_pre
                        buffer{i,j}(:,:,ind) = M_temp;
                    end
                end
            end
        end
        
%         ss = squeeze(shifts_up(:,:,t,:));
%         SS = mat2cell(ss,ones(size(ss,1),1),ones(size(ss,2),1),2);
%         Ic = mat2cell_ov(Yt,grid_size_fine,overlap_post);
%         M_fin = cellfun(@(x,s) shift_reconstruct(x,s,us_fac,1,1),Ic,SS,'un',0);

        
%         M_fin = cell(length(xx_fine)*length(yy_fine),1);
%         M_fin2 = M_fin;
%         parfor ii = 1:length(xx_fine)*length(yy_fine)
%             [i,j] = ind2sub([length(xx_fine),length(yy_fine)],ii)
%             true_grid = [xx_fine(i),min(xx_fine(i)+grid_size_fine(1)-1,d1),yy_fine(j),min(yy_fine(j)+grid_size_fine(2)-1,d2)];
%             extended_grid = [max(xx_fine(i)-overlap_post,1),min(xx_fine(i)+grid_size_fine(1)+overlap_post-1,d1),max(yy_fine(j)-overlap_post,1),min(yy_fine(j)+grid_size_fine(2)+overlap_post-1,d2)];
%             I_temp = Yt(extended_grid(1):extended_grid(2),extended_grid(3):extended_grid(4));
%             M_temp = shift_reconstruct(I_temp,shifts_up(i,j,t,:),us_fac,1,1);
%             %M_fin{ii} = M_temp(1+(true_grid(1)-extended_grid(1)):end-(extended_grid(2)-true_grid(2)),1+(true_grid(3)-extended_grid(3)):end-(extended_grid(4)-true_grid(4)));
%             M_fin2{ii} = M_temp;
%         end
%         ccc = 0;
%         for i = 1:length(yy_fine)
%             for j = 1:length(xx_fine)
%                 ccc = ccc + 1;
%                 if ~buffer_pre
%                     buffer{j,i}(:,:,ind) = M_fin2{ccc};
%                 end
%             end
%         end
%         M_fin2 = reshape(M_fin2,length(yy_fine),length(xx_fine));        
        %Mf = cell2mat(M_fin) - add_value;
        Mf = cell2mat_ov(M_fin,grid_size_fine,overlap_post,d1,d2) - add_value;
        
        if strcmpi(filetype,'mat');
            M_final(:,:,t) = Mf;
        elseif memmap
            M_final.M(:,:,t) = Mf;
        end                        
        if plot_flag && mod(t,2) == 0
            subplot(221); imagesc(Yt-add_value,[nnY,mmY]); title('Raw data','fontweight','bold','fontsize',14); 
                            xlabel(sprintf('Frame %i out of %i',t,T),'fontweight','bold','fontsize',14); set(gca,'Xtick',[],'Ytick',[]);
            subplot(222); imagesc(Mf,[nnY,mmY]);  title('Motion Corrected','fontweight','bold','fontsize',14); colormap('bone'); axis off;
            subplot(223); quiver(shifts_up(:,:,t,1),shifts_up(:,:,t,2),'Autoscale','off'); title('Motion vector field','fontweight','bold','fontsize',14); axis off;
            subplot(224); imagesc(cell2mat_ov(template,grid_size,overlap_pre,d1,d2)-add_value,[nnY,mmY]); title('Matching Template','fontweight','bold','fontsize',14); axis off
            drawnow;
            if make_avi  
                currFrame = getframe(fig);
                writeVideo(vidObj,currFrame);    
            end
        end
    
    end
if mot_uf == 1; shifts_up = shifts; end

if it == iter
    template = cellfun(@(x) x - add_value,template,'un',0);
    template = cell2mat_ov(template,grid_size,overlap_pre,d1,d2);
end
if ~strcmpi(filetype,'mat');
    M_final.shifts = shifts;
    M_final.shifts_up = shifts_up; 
    M_final.template = template;
end

if make_avi && plot_flag
    close(vidObj);
end
maxNumCompThreads('automatic');
end