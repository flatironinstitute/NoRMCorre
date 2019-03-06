function [M_final,shifts_g,template,options,col_shift] = normcorre_batch_even(Y,options,template)

% online motion correction through DFT subpixel registration
% Based on the dftregistration.m function from Manuel Guizar and Jim Fienup

% INPUTS
% Y:                Input data, can be already loaded in memory as a 3D
%                   tensor, a memory mapped file, or a pointer to a tiff stack
% options:          options structure for motion correction (optional, rigid registration is performed if not provided)
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
    if strcmpi(ext,'tif') || strcmpi(ext,'tiff')
        tiffInfo = imfinfo(Y);
        filetype = 'tif';
        T = length(tiffInfo);
        sizY = [tiffInfo(1).Height,tiffInfo(1).Width,T];
    elseif strcmpi(ext,'mat')
        filetype = 'mem';
        Y = matfile(Y,'Writable',true);
        details = whos(Y);
        var_sizes = [details.bytes];
        [~,var_ind] = max(var_sizes);
        var_name = details(var_ind).name;
        sizY = size(Y,var_name);
        T = sizY(end);
    elseif strcmpi(ext,'hdf5') || strcmpi(ext,'h5')
        filetype = 'hdf5';
        fileinfo = hdf5info(Y);
        data_name = fileinfo.GroupHierarchy.Datasets.Name;
        sizY = fileinfo.GroupHierarchy.Datasets.Dims;
        T = sizY(end);
    elseif strcmpi(ext,'raw')
        filetype = 'raw';
        fid = fopen(Y);
        FOV = [options.d1,options.d2];
        bitsize = options.bitsize;
        imsize = FOV(1)*FOV(2)*bitsize;                                                   % Bit size of single frame
        current_seek = ftell(fid);
        fseek(fid, 0, 1);
        file_length = ftell(fid);
        fseek(fid, current_seek, -1);
        T = file_length/imsize;
        sizY = [FOV,T];
        fclose(fid);        
    elseif strcmpi(ext,'avi')
        filetype = 'avi';
        sizY = size(read_file(Y));
        FOV = sizY(1:2);
        T = sizY(end);
    end    
elseif isobject(Y)
    filetype = 'mem';
    var_name = 'Y';
    sizY = size(Y,var_name);
    T = sizY(end);
else % array loaded in memory
    filetype = 'mat';
    sizY = size(Y);
    T = sizY(end);
end

nd = length(sizY)-1;                          % determine whether imaging is 2d or 3d
sizY = sizY(1:nd);
otherdims = repmat({':'},1,nd);
%% set default parameters if not present

if ~exist('options','var') || isempty(options)
    options = NoRMCorreSetParms('d1',sizY(1),'d2',sizY(2));
    if nd > 2; options.d3 = sizY(3); end
end

memmap = options.memmap;
grid_size = options.grid_size; 
mot_uf = options.mot_uf;
min_patch_size = options.min_patch_size;
min_diff = options.min_diff;
overlap_pre = options.overlap_pre;
overlap_post = options.overlap_post;
upd_template = options.upd_template;
bin_width = options.bin_width;
buffer_width = options.buffer_width;
max_dev_g = options.max_dev;
init_batch = options.init_batch;
us_fac = options.us_fac;
method = options.method;
filename = options.mem_filename;
iter = options.iter;
add_value = options.add_value;
max_shift = options.max_shift;
if strcmpi(options.boundary,'nan')
    fill_value = NaN;
else
    fill_value = add_value;
end
print_msg = options.print_msg;


while mod(T,bin_width) == 1
    if T == 1
        error('Movie appears to have only one frame. Use the function normcorre instead')        
    end
    bin_width = bin_width + 1;
end

%% first check for offset due to bi-directional scanning

if options.correct_bidir && isempty(options.col_shift)
    col_shift = correct_bidirectional_offset(Y,options.nFrames,options.bidir_us);
elseif ~isempty(options.col_shift)
    col_shift = options.col_shift;
else
    col_shift = 0;
end 
options.col_shift = col_shift;
if col_shift 
    if print_msg; fprintf('Offset %1.1d pixels due to bidirectional scanning detected. \n',col_shift); end
    if strcmpi(options.shifts_method,'fft')
        options.shifts_method = 'cubic';
        if print_msg; fprintf('Cubic shifts will be applied. \n'); end
    end
end
%% read initial batch and compute template

init_batch = min(T,init_batch);
interval = ceil(T/2-init_batch/2+1):floor(T/2+init_batch/2);
switch filetype
    case 'tif'
        Y_temp = read_file(Y,interval(1),init_batch,[],tiffInfo);
    case 'hdf5'
        Y_temp = read_file(Y,interval(1),init_batch);        
    case 'avi'
        Y_temp = read_file(Y,interval(1),init_batch);
    case 'mem'
        Y_temp = Y.(var_name)(otherdims{:},interval);
    case 'mat'
        Y_temp = Y(otherdims{:},interval);
    case 'raw'
        Y_temp = read_raw_file(Y,interval(1),init_batch,FOV,bitsize);
end
data_type = class(Y_temp);
Y_temp = single(Y_temp);
use_proj = true;
if nargin < 3 || isempty(template)
    if print_msg; fprintf('Registering the first %i frames just to obtain a good template....',init_batch); end
    template_in = median(Y_temp,nd+1)+add_value;
    fftTemp = fftn(template_in);
    for t = 1:size(Y_temp,nd+1)        
        if nd == 2
            [~,Greg] = dftregistration_min_max(fftTemp,fftn(Y_temp(:,:,t)),us_fac,-max_shift,max_shift,options.phase_flag);
        end
        if nd == 3 
            [~,Greg] = dftregistration_min_max_3d(fftTemp,fftn(Y_temp(:,:,:,t)),us_fac,-max_shift,max_shift,options.phase_flag); 
        end
        M_temp = real(ifftn(Greg));
        template_in = template_in*(t-1)/t + M_temp/t;
    end
    template_in = template_in + add_value;
    if print_msg; fprintf('..done. \n'); end
else
    template_in = single(template + add_value);
end

[d1,d2,d3,~] = size(Y_temp);
if nd == 2; d3 = 1; end
%% setup grids for patches
dim = [d1,d2,d3];
patches = construct_grid_even(grid_size,overlap_pre,dim,min_diff);
shifts_g = struct('shifts',cell(T,1),'shifts_up',cell(T,1),'diff',cell(T,1));
template_patches = split_frame(template_in,patches);
patch_size = size(template_patches);
n_patches = [length(unique(patches(:,1))),length(unique(patches(:,3))),length(unique(patches(:,5)))];
%%
%maxNumCompThreads(1);
temp_mat = template_in;

fftTemp =fft(fft(fft(template_patches,[],1),[],2),[],3);
fftTempMat = fftn(temp_mat);

if ~strcmpi(options.output_type,'mat')
    options.mem_batch_size = max(min(round(options.mem_batch_size/bin_width)*bin_width,T),1);
    mem_buffer = squeeze(zeros([dim,options.mem_batch_size],'single'));
end

switch lower(options.output_type)
    case 'mat'
        M_final = zeros([sizY,T],data_type);
    case 'memmap'
        M_final = matfile(filename,'Writable',true);
        if nd == 2; M_final.Y(d1,d2,T) = zeros(1,data_type); end
        if nd == 3; M_final.Y(d1,d2,d3,T) = zeros(1,data_type); end
        M_final.Yr(d1*d2*d3,T) = zeros(1,data_type);        
    case {'hdf5','h5'}
         if exist(options.h5_filename,'file')
            [pathstr,fname,ext] = fileparts(options.h5_filename);             
            new_filename = fullfile(pathstr,[fname,'_',datestr(now,30),ext]);
            warning_msg = ['File ',options.h5_filename,'already exists. Saving motion corrected file as',new_filename];            
            warning('%s',warning_msg);
            options.h5_filename = new_filename;
        end       
        M_final = options.h5_filename;
        if nd == 2
            h5create(options.h5_filename,['/',options.h5_groupname],[d1,d2,Inf],'Chunksize',[d1,d2,options.mem_batch_size],'Datatype',data_type);
        elseif nd == 3
            h5create(options.h5_filename,['/',options.h5_groupname],[d1,d2,d3,Inf],'Chunksize',[d1,d2,d3,options.mem_batch_size],'Datatype',data_type);
        end
    case {'tif','tiff'}
        M_final = ['motion corrected file has been saved as ', options.tiff_filename];
        opts_tiff.append = true;
        opts_tiff.big = true;
        if nd == 3
            error('Saving volumetric tiff stacks is currently not supported. Use a different filetype');
        end
    otherwise
        error('This filetype is currently not supported')
end   

cnt_buf = 0;
if print_msg; fprintf('Template initialization complete.  Now registering all the frames with new template. \n'); end
%%

prevstr = [];
for it = 1:iter
    for t = 1:bin_width:T
        switch filetype
            case 'tif'
                Ytm = single(read_file(Y, t, min(t+bin_width-1,T)-t+1, [], tiffInfo));
            case 'avi'
                Ytm = single(read_file(Y, t, min(t+bin_width-1,T)-t+1));
            case 'hdf5'
                Ytm = single(h5read(Y,data_name,[ones(1,nd),t],[sizY(1:nd),min(t+bin_width-1,T)-t+1]));
            case 'mem'
                Ytm = single(Y.(var_name)(otherdims{:},t:min(t+bin_width-1,T)));
            case 'mat'
                Ytm = single(Y(otherdims{:},t:min(t+bin_width-1,T)));
            case 'raw'
                Ytm = single(read_raw_file(Y,t,min(t+bin_width-1,T)-t+1,FOV,bitsize));                
        end
        
        Mf = zeros(size(Ytm),data_type);
        lY = size(Ytm,nd+1);
        shifts = struct('shifts',cell(lY,1),'shifts_up',cell(lY,1),'diff',cell(lY,1));
        for i = lY:-1:1
            Yt = Ytm(otherdims{:},i);
            future_results(i) = parfeval(@register_frame, 2, Yt,fftTempMat,fftTemp,patches,options);
        end
        for i = 1:lY
            [idx, shifts_temp, Mf_temp] = fetchNext(future_results);
            shifts(idx).shifts = shifts_temp;
            shifts(idx).shifts_up = shifts_temp;
            Mf(otherdims{:},idx) = Mf_temp;
        end
        
        shifts_g(t:min(t+bin_width-1,T)) = shifts;
        
        if ~strcmpi(options.output_type,'mat')
            rem_mem = rem(t+lY-1,options.mem_batch_size);
            if rem_mem == 0; rem_mem = options.mem_batch_size; end            
            mem_buffer(otherdims{:},rem_mem-lY+1:rem_mem) = cast(Mf,data_type);
        end
        if it == iter
            switch lower(options.output_type)
                case 'mat'
                    M_final(otherdims{:},t:min(t+bin_width-1,T)) = cast(Mf,data_type);
                case 'memmap'
                    if rem_mem == options.mem_batch_size || t+lY-1 == T
                        M_final.Y(otherdims{:},t+lY-rem_mem:t+lY-1) = mem_buffer(:,:,1:rem_mem);
                        M_final.Yr(:,t+lY-rem_mem:t+lY-1) = reshape(mem_buffer(1:d1*d2*d3*rem_mem),d1*d2*d3,rem_mem);
                    end      
                case {'hdf5','h5'}
                    if rem_mem == options.mem_batch_size || t+lY-1 == T
                        if nd == 2; h5write(options.h5_filename,['/',options.h5_groupname],mem_buffer(:,:,1:rem_mem),[ones(1,nd),t+lY-rem_mem],[sizY(1:nd),rem_mem]); end
                        if nd == 3; h5write(options.h5_filename,['/',options.h5_groupname],mem_buffer(:,:,:,1:rem_mem),[ones(1,nd),t+lY-rem_mem],[sizY(1:nd),rem_mem]); end
                    end
                case {'tif','tiff'}
                    if rem_mem == options.mem_batch_size || t+lY-1 == T
                        saveastiff(cast(mem_buffer(:,:,1:rem_mem),data_type),options.tiff_filename,opts_tiff);
                    end
            end        
        end
        if print_msg
            str = sprintf('%d out of %d frames registered, iteration %d out of %d..', t+lY-1, T, it, iter);
            refreshdisp(str, prevstr, t);
            prevstr=str; 
        end
        % update template
        if upd_template
            cnt_buf = cnt_buf + 1;                 
            if strcmpi(method{2},'mean')
                new_temp = nanmean(Mf,nd+1);
            elseif strcmpi(method{2},'median')
                new_temp = nanmedian(x,nd+1);
            end
            if any(isnan(new_temp(:)))            
                new_temp(isnan(new_temp)) =  template(isnan(new_temp));
            end
            if strcmpi(method{1},'mean')
                cnt = t/bin_width + 1;
                temp_mat = temp_mat*(cnt-1)/cnt + new_temp/cnt;
            elseif strcmpi(method{1},'median')
                ind_buffer = mod(cnt_buf,buffer_width) + buffer_width*(~mod(cnt_buf,buffer_width));
                buffer_med(otherdims{:},ind_buffer) = new_temp;
                temp_mat = nanmedian(buffer_med,nd+1);
            end
            template_patches = split_frame(temp_mat,patches);
            fftTemp = fft(fft(fft(template_patches,[],1),[],2),[],3);            
            fftTempMat = fftn(temp_mat);
        end
    end

    if it == iter
        template = temp_mat - add_value;
    end
    if memmap
        M_final.shifts = shifts_g;
        M_final.template = template;
    end
end
if print_msg; fprintf('\n'); end
maxNumCompThreads('automatic');