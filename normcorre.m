function [M_final,shifts,template,options] = normcorre(Y,options,template)

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

%% first determine filetype

if isa(Y,'char')
    [~,~,ext] = fileparts(Y);
    ext = ext(2:end);
    if strcmpi(ext,'tif') || strcmpi(ext,'tiff');
        tiffInfo = imfinfo(Y);
        filetype = 'tif';
        T = length(tiffInfo);
        sizY = [tiffInfo(1).Height,tiffInfo(1).Width,T];
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
    sizY = size(Y,'Y');
    T = sizY(end);
else % array loaded in memory
    filetype = 'mat';
    Y = single(Y);
    sizY = size(Y);
    T = sizY(end);
end

nd = length(sizY)-1;                          % determine whether imaging is 2d or 3d
sizY = sizY(1:nd);
%% set default parameters if not present

defoptions.memmap = false;                     % save motion corrected file in a mat file
if nd == 2
    defoptions.grid_size = [128,128,1];       % size of each patch to be corrected separately
elseif nd == 3
    defoptions.grid_size = [64,64,16];
end
if nd == 2
    defoptions.mot_uf = [4,4,1];              % upsampling factor within each patch
elseif nd == 3
    defoptions.mot_uf = [2,2,1];
end
defoptions.min_patch_size = [32,32,16];             % minimum patch size 
defoptions.overlap_pre = [16,16,2];                 % overlap between subsets within each patch
defoptions.overlap_post = [8,8,2];                  % overlap between subsets within each patch
defoptions.upd_template = true;                     % flag for updating template
defoptions.bin_width = 10;                          % width of buffer for computing the moving template
defoptions.buffer_width = 50;                       % number of local means to keep in memory
defoptions.init_batch = 30;                         % length of initial batch
defoptions.max_dev = [3,3,1];                       % maximum deviation around rigid translation
defoptions.us_fac = 5;                              % upsampling factor for subpixel registration
defoptions.method = {'median';'mean'};              % method for averaging the template
defoptions.plot_flag = false;                       % flag for plotting results in real time
defoptions.mem_filename = 'motion_corrected.mat';   % filename for motion corrected mat file
defoptions.use_parallel = false;                    % use parfor when breaking each frame into patches
defoptions.make_avi = false;                        % flag for making movie
defoptions.name = 'motion_corrected.avi';           % name of saved movie
defoptions.fr = 30;                                 % frame rate for saved movie
defoptions.iter = 1;                                % number of passes over the data
defoptions.add_value = 0;                           % add value to make dataset non-negative
%defoptions.write_tiff = false;               % save output as a tiff stack
%defoptions.out_name = 'motion_corrected.tif'; % name for output file name

if nargin == 1 || isempty(options); options = defoptions; end

if ~isfield(options,'memmap'); options.memmap = defoptions.memmap; end; memmap = options.memmap;
if ~isfield(options,'grid_size'); options.grid_size = defoptions.grid_size; end; grid_size = options.grid_size; 
if ~isfield(options,'mot_uf'); options.mot_uf = defoptions.mot_uf; end; mot_uf = options.mot_uf;
if ~isfield(options,'min_patch_size'); options.min_patch_size = defoptions.min_patch_size; end; min_patch_size = options.min_patch_size;
if ~isfield(options,'overlap_pre'); options.overlap_pre = defoptions.overlap_pre; end; overlap_pre = options.overlap_pre;
if ~isfield(options,'overlap_post'); options.overlap_post = defoptions.overlap_post; end; overlap_post = options.overlap_post;
if ~isfield(options,'upd_template'); options.upd_template = defoptions.upd_template; end; upd_template = options.upd_template;
if ~isfield(options,'bin_width'); options.bin_width = defoptions.bin_width; end; bin_width = options.bin_width;
if ~isfield(options,'buffer_width'); options.buffer_width = defoptions.buffer_width; end; buffer_width = options.buffer_width;
if ~isfield(options,'max_dev'); options.max_dev = defoptions.max_dev; end; max_dev = options.max_dev;
if ~isfield(options,'init_batch'); options.init_batch = defoptions.init_batch; end; init_batch = options.init_batch;
if ~isfield(options,'us_fac'); options.us_fac = defoptions.us_fac; end; us_fac = options.us_fac;
if ~isfield(options,'method'); options.method = defoptions.method; end; method = options.method;
if ~isfield(options,'plot_flag'); options.plot_flag = defoptions.plot_flag; end; plot_flag = options.plot_flag*(nd==2);
if ~isfield(options,'mem_filename'); options.mem_filename = defoptions.mem_filename; end; filename = options.mem_filename;
if ~isfield(options,'use_parallel'); options.use_parallel = defoptions.use_parallel; end; use_parallel = options.use_parallel;
if ~isfield(options,'make_avi'); options.make_avi = defoptions.make_avi; end; make_avi = options.make_avi;
if ~isfield(options,'name'); options.name = defoptions.name; end; name = options.name;
if ~isfield(options,'fr'); options.fr = defoptions.fr; end; fr = options.fr;
if ~isfield(options,'iter'); options.iter = defoptions.iter; end; iter = options.iter;
if ~isfield(options,'add_value'); options.add_value = defoptions.add_value; end; add_value = options.add_value;
if isscalar(grid_size); grid_size = grid_size*ones(1,nd); end; if length(grid_size) == 2; grid_size(3) = 1; end
if isscalar(mot_uf); mot_uf = mot_uf*ones(1,nd); end; if length(mot_uf) == 2; mot_uf(3) = 1; end
if isscalar(overlap_pre); overlap_pre = overlap_pre*ones(1,nd); end; if length(overlap_pre) == 2; overlap_pre(3) = 1; end
if isscalar(overlap_post); overlap_post = overlap_post*ones(1,nd); end; if length(overlap_post) == 2; overlap_post(3) = 1; end
if isscalar(max_dev); max_dev = max_dev*ones(1,nd); end; if length(max_dev) == 2; max_dev(3) = 1; end

if ~isfield(options,'max_shift'); 
    max_shift = grid_size./mot_uf; 
else
    if isscalar(options.max_shift)
        options.max_shift = options.max_shift*ones(1,3);
    end
    max_shift = min(options.max_shift,grid_size./mot_uf);
end

%% read initial batch and compute template
init_batch = min(T,init_batch);
perm = randperm(T,init_batch);
switch filetype
    case 'tif'
        Y1 = imread(Y,'Index',1,'Info',tiffInfo);
        Y_temp = zeros(sizY(1),sizY(2),init_batch,'like',Y1);
        Y_temp(:,:,1) = Y1;
        for tt = 1:init_batch
            Y_temp(:,:,tt) = imread(Y,'Index',perm(tt),'Info',tiffInfo);
        end
    case 'hdf5'
        Y_temp = bigread2(Y,1,init_batch);        
    case 'mem'
        if nd == 2; Y_temp = Y.Y(:,:,1:init_batch); elseif nd == 3; Y_temp = Y.Y(:,:,:,1:init_batch); end
    case 'mat'
        if nd == 2; Y_temp = Y(:,:,perm); elseif nd == 3; Y_temp = Y(:,:,:,perm); end
end

data_type = class(Y_temp);
Y_temp = single(Y_temp);

if nargin < 3 || isempty(template)
    template_in = median(Y_temp,nd+1)+add_value;
else
    template_in = template + add_value;
end

[d1,d2,d3,~] = size(Y_temp);
if nd == 2; d3 = 1; end
%% setup grids for patches

[xx_s,xx_f,yy_s,yy_f,zz_s,zz_f,xx_us,xx_uf,yy_us,yy_uf,zz_us,zz_uf] = construct_grid(grid_size,mot_uf,d1,d2,d3,min_patch_size);
shifts = struct('shifts',cell(T,1),'shifts_up',cell(T,1));
temp_cell = mat2cell_ov(template_in,xx_us,xx_uf,yy_us,yy_uf,zz_us,zz_uf,overlap_post,sizY);

%% precompute some quantities that are used repetitively for template matching and applying shifts
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
            extended_grid = [max(xx_us(i)-overlap_post(1),1),min(xx_uf(i)+overlap_post(1),d1),max(yy_us(j)-overlap_post(2),1),min(yy_uf(j)+overlap_post(2),d2),max(zz_us(k)-overlap_post(3),1),min(zz_uf(k)+overlap_post(3),d3)];            
            Bs{i,j,k} = permute(construct_weights([xx_us(i),xx_uf(i),yy_us(j),yy_uf(j),zz_us(k),zz_uf(k)],extended_grid),[2,1,3]); 
        end
    end
end
if nd == 2; Np = cellfun(@(x) 0,Nr,'un',0); end

%%
%maxNumCompThreads(2);
template = mat2cell_ov(template_in,xx_s,xx_f,yy_s,yy_f,zz_s,zz_f,overlap_pre,sizY);
temp_mat = template_in;
fftTemp = cellfun(@fftn,template,'un',0);
fftTempMat = fftn(temp_mat);
if nd == 2; buffer = mat2cell_ov(zeros(d1,d2,bin_width,'single'),xx_s,xx_f,yy_s,yy_f,zz_s,zz_f,overlap_pre,sizY); end
if nd == 3; buffer = mat2cell_ov(zeros(d1,d2,d3,bin_width,'single'),xx_s,xx_f,yy_s,yy_f,zz_s,zz_f,overlap_pre,sizY); end


if ~strcmpi(options.output_type,'mat')
    %options.mem_batch_size = min(round(options.mem_batch_size/bin_width)*bin_width,T);
    if nd == 2; mem_buffer = zeros(d1,d2,options.mem_batch_size,'single'); end
    if nd == 3; mem_buffer = zeros(d1,d2,d3,options.mem_batch_size,'single'); end
end

switch lower(options.output_type)
    case 'mat'
        M_final = zeros([sizY,T]);
    case 'memmap'
        M_final = matfile(filename,'Writable',true);
        if nd == 2; M_final.Y(d1,d2,T) = zeros(1,data_type); end
        if nd == 3; M_final.Y(d1,d2,d3,T) = zeros(1,data_type); end
        M_final.Yr(d1*d2*d3,T) = zeros(1,data_type);        
    case {'hdf5','h5'}
        M_final = ['motion corrected file has been saved as ', options.h5_filename];
        if nd == 2
            h5create(options.h5_filename,['/',options.h5_groupname],[d1,d2,Inf],'Chunksize',[d1,d2,options.mem_batch_size],'Datatype',data_type);
        elseif nd == 3
            h5create(options.h5_filename,['/',options.h5_groupname],[d1,d2,d3,Inf],'Chunksize',[d1,d2,d3,options.mem_batch_size],'Datatype',data_type);
        end
    otherwise
        error('This filetype is currently not supported')
end   

%%
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
%%
for it = 1:iter
    if it < iter; plot_flag = 0; else plot_flag = options.plot_flag; end
    for t = 1:T
        switch filetype
            case 'tif'
                Yt = single(imread(Y,'Index',t,'Info',tiffInfo));
            case 'hdf5'
                Yt = single(h5read(Y,'/mov',[ones(1,nd),t],[sizY(1:nd),1]));
            case 'mem'
                if nd == 2; Yt = single(Y.Y(:,:,t)); end
                if nd == 3; Yt = single(Y.Y(:,:,:,t)); end
            case 'mat'
                if nd == 2; Yt = single(Y(:,:,t)); end
                if nd == 3; Yt = single(Y(:,:,:,t)); end
        end        
        minY = min(Yt(:));
        maxY = max(Yt(:));
        Yt = Yt + add_value;
        ind = rem(t,bin_width) + bin_width*(rem(t,bin_width)==0);
        Yc = mat2cell_ov(Yt,xx_s,xx_f,yy_s,yy_f,zz_s,zz_f,overlap_pre,sizY);
        fftY = cellfun(@fftn, Yc, 'un',0);
        
        M_fin = cell(length(xx_us),length(yy_us),length(zz_us)); %zeros(size(Y_temp));
        shifts_temp = zeros(length(xx_s),length(yy_s),length(zz_s),nd); 
        diff_temp = zeros(length(xx_s),length(yy_s),length(zz_s));
        if numel(M_fin) > 1           
            if nd == 2; out_rig = dftregistration_min_max(fftTempMat,fftn(Yt),us_fac,-max_shift,max_shift,options.phase_flag); lb = out_rig(3:4); ub = out_rig(3:4); end
            if nd == 3; out_rig = dftregistration_min_max_3d(fftTempMat,fftn(Yt),us_fac,-max_shift,max_shift,options.phase_flag); lb = out_rig(3:5); ub = out_rig(3:5); end
        else
            lb = -max_shift(1,nd);
            ub = max_shift(1,nd);
            max_dev = 0*max_dev;
        end
        if ~use_parallel
            for i = 1:length(xx_s)
                for j = 1:length(yy_s)           
                    for k = 1:length(zz_s)
                        if nd == 2
                            %[output,Greg] = dftregistration_max(fftTemp{i,j,k},fftY{i,j,k},us_fac,max_shift);        
                            [output,Greg] = dftregistration_min_max(fftTemp{i,j,k},fftY{i,j,k},us_fac,lb-max_dev(1:2),ub+max_dev(1:2),options.phase_flag);  
                        elseif nd == 3
                            %[output,Greg] = dftregistration_max_3d(fftTemp{i,j,k},fftY{i,j,k},us_fac,max_shift);
                            [output,Greg] = dftregistration_min_max_3d(fftTemp{i,j,k},fftY{i,j,k},us_fac,lb-max_dev,ub+max_dev,options.phase_flag); 
                            shifts_temp(i,j,k,3) = output(5);
                        end
                        M_temp = real(ifftn(Greg));
                        M_temp = remove_boundaries(M_temp,output(3:end),'copy',template{i,j,k});
                        if nd == 2; buffer{i,j,k}(:,:,ind) = M_temp;  end    
                        if nd == 3; buffer{i,j,k}(:,:,:,ind) = M_temp;  end 
                        shifts_temp(i,j,k,1) = output(3);
                        shifts_temp(i,j,k,2) = output(4); 
                        diff_temp(i,j,k) = output(2);
                        if mot_uf == 1
                            M_fin{i,j,k} = remove_boundaries(M_temp,output(3:end),options.boundary,template{i,j,k},add_value);
                        end                                               
                    end
                end
            end            
        else
            Mt2 = cell(length(xx_s)*length(yy_s)*length(zz_s),1);            
            shifts_cell = cell(length(xx_s)*length(yy_s)*length(zz_s),1); 
            diff_cell = cell(length(xx_s)*length(yy_s)*length(zz_s),1); 
            parfor ii = 1:length(xx_s)*length(yy_s)*length(zz_s)
                [i,j,k] = ind2sub([length(xx_s),length(yy_s),length(zz_s)],ii)
                %if nd == 2; [output,Greg] = dftregistration_max(fftTemp{i,j,k},fftY{i,j,k},us_fac,max_shift); end
                %if nd == 3; [output,Greg] = dftregistration_max_3d(fftTemp{i,j,k},fftY{i,j,k},us_fac,max_shift); end
                if nd == 2; [output,Greg] = dftregistration_min_max(fftTemp{i,j,k},fftY{i,j,k},us_fac,lb-max_dev(1:2),ub+max_dev(1:2),options.phase_flag); end
                if nd == 3; [output,Greg] = dftregistration_min_max_3d(fftTemp{i,j,k},fftY{i,j,k},us_fac,lb-max_dev,ub+max_dev,options.phase_flag); end                
                M_temp = real(ifftn(Greg));
                Mt2{ii} = M_temp;
                shifts_cell{ii} = output(3:end);
                diff_cell{ii} = output(2);
            end
            for ii = 1:length(xx_s)*length(yy_s)*length(zz_s)
                 [i,j,k] = ind2sub([length(xx_s),length(yy_s),length(zz_s)],ii);
                 buffer{i,j,k}(:,:,ind) = Mt2{ii};
                 if mot_uf == 1
                     M_fin{i,j,k} = Mt2{ii};
                 end
                 shifts_temp(i,j,k,:) = shifts_cell{ii};
                 diff_temp(i,j,k) = diff_cell{ii};
            end            
        end     
        shifts(t).shifts = squeeze(shifts_temp);
        shifts(t).diff = diff_temp;
        
        if any(mot_uf > 1)
            shifts_up = imresize(shifts_temp,[length(xx_uf),length(yy_uf)]);
            diff_up = imresize(diff_temp,[length(xx_uf),length(yy_uf)]);
%             shifts_up(:,:,:,1) = bound_shifts(shifts_up(:,:,:,1),0.5);
%             shifts_up(:,:,:,2) = bound_shifts(shifts_up(:,:,:,2),0.5);
            %if 
            if mot_uf(3) > 1
                shifts_up = reshape(imresize(reshape(shifts_up,[length(xx_uf)*length(yy_uf),length(zz_f),nd]),[length(xx_uf)*length(yy_uf),length(zz_uf)]),[length(xx_uf),length(yy_uf),length(zz_uf),nd]);
                diff_up = reshape(imresize(reshape(diff_up,[length(xx_uf)*length(yy_uf),length(zz_f)]),[length(xx_uf)*length(yy_uf),length(zz_uf)]),[length(xx_uf),length(yy_uf),length(zz_uf)]);
            end
            shifts(t).shifts_up = shifts_up;
            shifts(t).diff = diff_up;
            for i = 1:length(xx_uf)
                for j = 1:length(yy_uf)
                    for k = 1:length(zz_uf)
                        extended_grid = [max(xx_us(i)-overlap_post(1),1),min(xx_uf(i)+overlap_post(1),d1),max(yy_us(j)-overlap_post(2),1),min(yy_uf(j)+overlap_post(2),d2),max(zz_us(k)-overlap_post(3),1),min(zz_uf(k)+overlap_post(3),d3)];
                        I_temp = Yt(extended_grid(1):extended_grid(2),extended_grid(3):extended_grid(4),extended_grid(5):extended_grid(6));
                        M_fin{i,j,k} = shift_reconstruct(I_temp,shifts_up(i,j,k,:),diff_up(i,j,k),us_fac,Nr{i,j,k},Nc{i,j,k},Np{i,j,k},options.boundary,add_value);
                    end
                end
            end
        else
            shifts_up = shifts_temp;
            shifts(t).shifts_up = shifts(t).shifts;
        end
%         ss = squeeze(shifts_up(:,:,t,:));
%         SS = mat2cell(ss,ones(size(ss,1),1),ones(size(ss,2),1),2);
%         Ic = mat2cell_ov(Yt,grid_size_fine,overlap_post);

        gx = max(abs(reshape(diff(shifts_up,[],1),[],1)));
        gy = max(abs(reshape(diff(shifts_up,[],2),[],1)));
        gz = max(abs(reshape(diff(shifts_up,[],3),[],1)));
        flag_interp = max([gx;gy;gz;0])<0.5;      % detect possible smearing

        if flag_interp    
            Mf = cell2mat_ov_sum(M_fin,xx_us,xx_uf,yy_us,yy_uf,zz_us,zz_uf,overlap_post,sizY,Bs) - add_value;
        else            
            Mf = cell2mat_ov(M_fin,xx_us,xx_uf,yy_us,yy_uf,zz_us,zz_uf,overlap_post,sizY) - add_value;
        end
        
        Mf(Mf<minY) = minY;
        Mf(Mf>maxY) = maxY;
        if ~strcmpi(options.output_type,'mat')
            rem_mem = rem(t,options.mem_batch_size);
            if rem_mem == 0; rem_mem = options.mem_batch_size; end            
            if nd == 2; mem_buffer(:,:,rem_mem) = cast(Mf,data_type); end
            if nd == 3; mem_buffer(:,:,:,rem_mem) = cast(Mf,data_type); end
        end
        switch lower(options.output_type)
            case 'mat'
                if nd == 2; M_final(:,:,t) = Mf; end
                if nd == 3; M_final(:,:,:,t) = Mf; end
            case 'memmap'
                if rem_mem == options.mem_batch_size || t == T
                    if nd == 2; M_final.Y(:,:,t-rem_mem+1:t) = mem_buffer(:,:,1:rem_mem); end
                    if nd == 3; M_final.Y(:,:,:,t-rem_mem+1:t) = mem_buffer(:,:,:,1:rem_mem); end
                    M_final.Yr(:,t-rem_mem+1:t) = reshape(mem_buffer(1:d1*d2*d3*rem_mem),d1*d2*d3,rem_mem);
                end      
            case {'hdf5','h5'}
                if rem_mem == options.mem_batch_size || t== T
                    if nd == 2; h5write(options.h5_filename,['/',options.h5_groupname],mem_buffer(:,:,1:rem_mem),[ones(1,nd),t-rem_mem+1],[sizY(1:nd),rem_mem]); end
                    if nd == 3; h5write(options.h5_filename,['/',options.h5_groupname],mem_buffer(:,:,:,1:rem_mem),[ones(1,nd),t-rem_mem+1],[sizY(1:nd),rem_mem]); end
                end
        end         
        
        if mod(t,bin_width) == 0 && upd_template
            fprintf('%i out of %i frames registered, iteration %i out of %i \n',t,T,it,iter)
            cnt_buf = cnt_buf + 1;                
            if strcmpi(method{2},'mean')
                new_temp = cellfun(@(x) nanmean(x,nd+1), buffer, 'UniformOutput',false);
            elseif strcmpi(method{2},'median');
                new_temp = cellfun(@(x) nanmedian(x,nd+1), buffer, 'UniformOutput', false);
            end
            if strcmpi(method{1},'mean')
                cnt = t/bin_width + 1;
                template = cellfun(@plus, cellfun(@(x) x*(cnt-1)/cnt, template,'un',0), cellfun(@(x) x*1/cnt, new_temp,'un',0), 'un',0);
            elseif strcmpi(method{1},'median');
                if cnt_buf <= buffer_width
                    if nd == 2; buffer_med(:,:,cnt_buf) = cell2mat_ov(new_temp,xx_s,xx_f,yy_s,yy_f,zz_s,zz_f,overlap_pre,sizY); end
                    if nd == 3; buffer_med(:,:,:,cnt_buf) = cell2mat_ov(new_temp,xx_s,xx_f,yy_s,yy_f,zz_s,zz_f,overlap_pre,sizY); end
                else
                    buffer_med = circshift(buffer_med,[zeros(1,nd),-1]);
                    if nd == 2; buffer_med(:,:,buffer_width) = cell2mat_ov(new_temp,xx_s,xx_f,yy_s,yy_f,zz_s,zz_f,overlap_pre,sizY); end
                    if nd == 3; buffer_med(:,:,:,buffer_width) = cell2mat_ov(new_temp,xx_s,xx_f,yy_s,yy_f,zz_s,zz_f,overlap_pre,sizY); end
                end
                template = mat2cell_ov(nanmedian(buffer_med,nd+1),xx_s,xx_f,yy_s,yy_f,zz_s,zz_f,overlap_pre,sizY);
            end
            fftTemp = cellfun(@fftn, template, 'un',0);
            temp_mat = cell2mat_ov(template,xx_s,xx_f,yy_s,yy_f,zz_s,zz_f,overlap_pre,sizY);
            fftTempMat = fftn(temp_mat);
        end        
        
        if plot_flag && mod(t,1) == 0
            subplot(221); imagesc(Yt-add_value,[nnY,mmY]); title('Raw data','fontweight','bold','fontsize',14); 
                            xlabel(sprintf('Frame %i out of %i',t,T),'fontweight','bold','fontsize',14); set(gca,'Xtick',[],'Ytick',[]);
            subplot(222); imagesc(Mf,[nnY,mmY]);  title('Motion Corrected','fontweight','bold','fontsize',14); colormap('bone'); axis off;
            subplot(223); quiver(shifts_up(:,:,:,1),shifts_up(:,:,:,2),'Autoscale','off'); title('Motion vector field','fontweight','bold','fontsize',14); axis off;
            subplot(224); imagesc(cell2mat_ov(template,xx_s,xx_f,yy_s,yy_f,zz_s,zz_f,overlap_pre,sizY)-add_value,[nnY,mmY]); title('Matching Template','fontweight','bold','fontsize',14); axis off
            drawnow;
            if make_avi  
                currFrame = getframe(fig);
                writeVideo(vidObj,currFrame);    
            end
        end   
    end

if it == iter
    template = cellfun(@(x) x - add_value,template,'un',0);
    template = cell2mat_ov(template,xx_s,xx_f,yy_s,yy_f,zz_s,zz_f,overlap_pre,sizY);
end
if memmap
    M_final.shifts = shifts;
    M_final.template = template;
end

if make_avi && plot_flag
    close(vidObj);
end
maxNumCompThreads('automatic');
end