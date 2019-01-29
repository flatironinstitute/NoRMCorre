function [M_final,shifts,template,options,col_shift] = normcorre(Y,options,template)

% online motion correction through DFT subpixel registration
% Based on the dftregistration.m function from Manuel Guizar and Jim Fienup

% INPUTS
% Y:                Input data, can be already loaded in memory as a 3D
%                   tensor, a memory mapped file, or a pointer to a tiff stack
% options:          options structure for motion correction (optional, rigid registration is performed if not provided)
% template:         provide template (optional)

% OUTPUTS
% M_final:          motion corrected data
% shifts:           originally calculated shifts
% template:         calculated template
% options:          options structure (if modified)
% col_shift:        relative shift due to bi-directional scanning

%% first determine filetype

nd = 2 + (options.d3 > 1); %max(length(sizY)-1,2);                    % determine whether imaging is 2d or 3d

if isa(Y,'char')
    [~,~,ext] = fileparts(Y);
    ext = ext(2:end);
    if strcmpi(ext,'tif') || strcmpi(ext,'tiff')
        tiffInfo = imfinfo(Y);
        filetype = 'tif';
        T = length(tiffInfo);
        if nd == 3
            sizY = [tiffInfo(1).Height,tiffInfo(1).Width,T,1];
        else
            sizY = [tiffInfo(1).Height,tiffInfo(1).Width,T];
        end
    elseif strcmpi(ext,'mat')
        filetype = 'mem';
        Y = matfile(Y,'Writable',true);
        details = whos(Y);
        var_sizes = [details.bytes];
        [~,var_ind] = max(var_sizes);
        var_name = details(var_ind).name;
        sizY = size(Y,var_name);
    elseif strcmpi(ext,'hdf5') || strcmpi(ext,'h5')
        filetype = 'hdf5';
        fileinfo = hdf5info(Y);
        sizY = fileinfo.GroupHierarchy.Datasets.Dims;
    end    
elseif isobject(Y)
    filetype = 'mem';
    var_name = 'Y';
    sizY = size(Y,var_name);
else % array loaded in memory
    filetype = 'mat';
    Y = single(Y);
    sizY = size(Y);    
end

if length(sizY) == nd
    T = 1;
else
    T = sizY(nd+1);
end
sizY = sizY(1:nd);

%% set default parameters if not present

if ~exist('options','var') || isempty(options)    
    options = NoRMCorreSetParms('d1',sizY(1),'d2',sizY(2));
    if nd > 2; options.d3 = sizY(3); end
end

options.bin_width = min(options.bin_width,T+1);
options.mem_batch_size = min(options.mem_batch_size,T);
options.buffer_width = min(options.buffer_width,ceil(T/options.bin_width));
options.init_batch = min(options.init_batch,T);

memmap = options.memmap;
grid_size = options.grid_size; 
mot_uf = options.mot_uf;
min_patch_size = options.min_patch_size;
overlap_pre = options.overlap_pre;
overlap_post = options.overlap_post;
upd_template = options.upd_template;
bin_width = options.bin_width;
buffer_width = options.buffer_width;
max_dev = options.max_dev;
init_batch = options.init_batch;
us_fac = options.us_fac;
method = options.method;
plot_flag = options.plot_flag*(nd==2);
filename = options.mem_filename;
use_parallel = options.use_parallel;
make_avi = options.make_avi;
name = options.name;
fr = options.fr;
iter = options.iter;
add_value = options.add_value;
max_shift = options.max_shift;
print_msg = options.print_msg;
if strcmpi(options.boundary,'nan')
    fill_value = NaN;
else
    fill_value = add_value;
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
if exist('template','var')
    init_batch = min(init_batch,1);
end
switch filetype
    case 'tif'
        Y_temp = read_file(Y,interval(1),init_batch,[],tiffInfo);
    case 'hdf5'
        Y_temp = read_file(Y,interval(1),init_batch);        
    case 'mem'
        if nd == 2; Y_temp = Y.(var_name)(:,:,interval); elseif nd == 3; Y_temp = Y.(var_name)(:,:,:,interval); end
    case 'mat'
        if nd == 2; Y_temp = Y(:,:,interval); elseif nd == 3; Y_temp = Y(:,:,:,interval); end
end

data_type = class(Y_temp);
Y_temp = single(Y_temp);

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

[xx_s,xx_f,yy_s,yy_f,zz_s,zz_f,xx_us,xx_uf,yy_us,yy_uf,zz_us,zz_uf] = construct_grid(grid_size,mot_uf,d1,d2,d3,min_patch_size);
shifts = struct('shifts',cell(T,1),'shifts_up',cell(T,1),'diff',cell(T,1));
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
        M_final = options.tiff_filename;
        opts_tiff.append = true;
        opts_tiff.big = true;
        if nd == 3
            error('Saving volumetric tiff stacks is currently not supported. Use a different filetype');
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
if print_msg; fprintf('Template initialization complete.  Now registering all the frames with new template. \n'); end
%%
prevstr = [];
for it = 1:iter
    if it < iter; plot_flag = 0; else plot_flag = options.plot_flag; end
    for t = 1:T
        switch filetype
            case 'tif'
                Yt = single(imread(Y,'Index',t,'Info',tiffInfo));
            case 'hdf5'
                Yt = single(h5read(Y,'/mov',[ones(1,nd),t],[sizY(1:nd),1]));
            case 'mem'
                if nd == 2; Yt = single(Y.(var_name)(:,:,t)); end
                if nd == 3; Yt = single(Y.(var_name)(:,:,:,t)); end
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
            if nd == 3; out_rig = dftregistration_min_max_3d(fftTempMat,fftn(Yt),1,-max_shift,max_shift,options.phase_flag); lb = out_rig(3:5); ub = out_rig(3:5); end
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
                        if all([length(xx_s),length(yy_s),length(zz_s)] == 1)
                            M_fin{i,j,k} = remove_boundaries(M_temp,output(3:end),options.boundary,template{i,j,k},add_value);
                        end                                               
                    end
                end
            end      
        else
            Mt2 = cell(length(xx_s)*length(yy_s)*length(zz_s),1);            
            shifts_cell = cell(length(xx_s)*length(yy_s)*length(zz_s),1); 
            diff_cell = cell(length(xx_s)*length(yy_s)*length(zz_s),1); 
            for ii = length(xx_s)*length(yy_s)*length(zz_s):-1:1
                [i,j,k] = ind2sub([length(xx_s),length(yy_s),length(zz_s)],ii);
                if nd == 2; future_results(ii) = parfeval(@dftregistration_min_max,2,fftTemp{i,j,k},fftY{i,j,k},us_fac,lb-max_dev(1:2),ub+max_dev(1:2),options.phase_flag); end
                if nd == 3; future_results(ii) = parfeval(@dftregistration_min_max_3d,2,fftTemp{i,j,k},fftY{i,j,k},us_fac,lb-max_dev,ub+max_dev,options.phase_flag); end
            end
            for i = 1:length(xx_s)*length(yy_s)*length(zz_s)
                [ii,output,Greg] = fetchNext(future_results);
                M_temp = real(ifftn(Greg));
                Mt2{ii} = M_temp;
                shifts_cell{ii} = output(3:end);
                diff_cell{ii} = output(2);
            end
%             parfor ii = 1:length(xx_s)*length(yy_s)*length(zz_s)
%                 [i,j,k] = ind2sub([length(xx_s),length(yy_s),length(zz_s)],ii);
%                 %if nd == 2; [output,Greg] = dftregistration_max(fftTemp{i,j,k},fftY{i,j,k},us_fac,max_shift); end
%                 %if nd == 3; [output,Greg] = dftregistration_max_3d(fftTemp{i,j,k},fftY{i,j,k},us_fac,max_shift); end
%                 if nd == 2; [output,Greg] = dftregistration_min_max(fftTemp{i,j,k},fftY{i,j,k},us_fac,lb-max_dev(1:2),ub+max_dev(1:2),options.phase_flag); end
%                 if nd == 3; [output,Greg] = dftregistration_min_max_3d(fftTemp{i,j,k},fftY{i,j,k},us_fac,lb-max_dev,ub+max_dev,options.phase_flag); end                
%                 M_temp = real(ifftn(Greg));
%                 Mt2{ii} = M_temp;
%                 shifts_cell{ii} = output(3:end);
%                 diff_cell{ii} = output(2);
%             end
            for ii = 1:length(xx_s)*length(yy_s)*length(zz_s)
                 [i,j,k] = ind2sub([length(xx_s),length(yy_s),length(zz_s)],ii);
                 if nd == 2; buffer{i,j,k}(:,:,ind) = Mt2{ii}; end
                 if nd == 3; buffer{i,j,k}(:,:,:,ind) = Mt2{ii}; end
                 if mot_uf == 1
                     M_fin{i,j,k} = Mt2{ii};
                 end
                 shifts_temp(i,j,k,:) = shifts_cell{ii};
                 diff_temp(i,j,k) = diff_cell{ii};
            end            
        end     
        shifts(t).shifts = shifts_temp;
        shifts(t).diff = diff_temp;
        
        switch lower(options.shifts_method)
            case 'fft'
                if any([length(xx_s),length(yy_s),length(zz_s)] > 1)
                    if ~isfield(options,'shifts_method'); options.shifts_method = 'FFT'; end                                         
                    if mot_uf(3) > 1                
                        do = [length(xx_us),length(yy_us),length(zz_us)]./[length(xx_s),length(yy_s),length(zz_s)];
                        ds = [length(xx_s),length(yy_s),length(zz_s)];
                        dim = [length(xx_us),length(yy_us),length(zz_us)];
                        [Xq,Yq,Zq] = meshgrid(linspace((1+1/do(2))/2,ds(2)+(1-1/do(2))/2,dim(2)),linspace((1+1/do(1))/2,ds(1)+(1-1/do(1))/2,dim(1)),linspace((1+1/do(3))/2,ds(3)+(1-1/do(3))/2,dim(3)));
                        %tform = affine3d(diag([mot_uf(:);1]));
                        %tform = affine3d(diag([mot_uf([2,1,3])';1]));
                        %diff_up = imwarp(diff_temp,tform,'OutputView',imref3d([length(xx_uf),length(yy_uf),length(zz_uf)]));
                        %diff_up = imwarp(diff_temp,tform,'OutputView',imref3d([length(xx_uf),length(yy_uf),length(zz_uf)]),'SmoothEdges',true);
                        diff_up = interp3(diff_temp,Xq,Yq,Zq,'makima');
                        shifts_up = zeros([size(diff_up),3]);
                        %for dm = 1:3; shifts_up(:,:,:,dm) = imwarp(shifts_temp(:,:,:,dm),tform,'OutputView',imref3d([length(xx_uf),length(yy_uf),length(zz_uf)]),'SmoothEdges',true); end
                        for dm = 1:3; shifts_up(:,:,:,dm) = interp3(shifts_temp(:,:,:,dm),Xq,Yq,Zq,'makima'); end
                    else
                        shifts_up = imresize(shifts_temp,[length(xx_uf),length(yy_uf)]);
                        diff_up = imresize(diff_temp,[length(xx_uf),length(yy_uf)]);
                    end

                    shifts(t).shifts_up = shifts_up;
                    shifts(t).diff = diff_up;
                    for i = 1:length(xx_uf)
                        for j = 1:length(yy_uf)
                            for k = 1:length(zz_uf)
                                extended_grid = [max(xx_us(i)-overlap_post(1),1),min(xx_uf(i)+overlap_post(1),d1),max(yy_us(j)-overlap_post(2),1),min(yy_uf(j)+overlap_post(2),d2),max(zz_us(k)-overlap_post(3),1),min(zz_uf(k)+overlap_post(3),d3)];
                                I_temp = Yt(extended_grid(1):extended_grid(2),extended_grid(3):extended_grid(4),extended_grid(5):extended_grid(6));
                                M_fin{i,j,k} = shift_reconstruct(I_temp,shifts_up(i,j,k,:),diff_up(i,j,k),us_fac,Nr{i,j,k},Nc{i,j,k},Np{i,j,k},options.boundary,add_value);                        
                                %M_fin{i,j,k} = shift_reconstruct2(I_temp,shifts_up(i,j,k,:),'bilinear',diff_up(i,j,k),us_fac,Nr{i,j,k},Nc{i,j,k},Np{i,j,k},options.boundary,add_value);
                            end
                        end
                    end                        
                else
                    shifts_up = shifts_temp;
                    shifts(t).shifts_up = shifts(t).shifts;
                end
                gx = max(abs(reshape(diff(shifts_up,[],1),[],1)));
                gy = max(abs(reshape(diff(shifts_up,[],2),[],1)));
                gz = max(abs(reshape(diff(shifts_up,[],3),[],1)));
                flag_interp = max([gx;gy;gz;0])<0.5;      % detect possible smearing

                if flag_interp    
                    Mf = cell2mat_ov_sum(M_fin,xx_us,xx_uf,yy_us,yy_uf,zz_us,zz_uf,overlap_post,sizY,Bs) - add_value;
                else            
                    Mf = cell2mat_ov(M_fin,xx_us,xx_uf,yy_us,yy_uf,zz_us,zz_uf,overlap_post,sizY) - add_value;
                end                             
        
            otherwise
                shifts(t).shifts_up = shifts(t).shifts;
                if nd == 3                
                    shifts_up = zeros([options.d1,options.d2,options.d3,3]);
                    do = size(shifts_up)./size(shifts_temp);
                    ds = size(shifts_temp);
                    dim = [options.d1,options.d2,options.d3];
                    if (0)
                        [Xq,Yq,Zq] = meshgrid(linspace(1,ds(2),dim(2)),linspace(1,ds(1),dim(1)),linspace(1,ds(3),dim(3)));
                    else
                        [Xq,Yq,Zq] = meshgrid(linspace((1+1/do(2))/2,ds(2)+(1-1/do(2))/2,dim(2)),linspace((1+1/do(1))/2,ds(1)+(1-1/do(1))/2,dim(1)),linspace((1+1/do(3))/2,ds(3)+(1-1/do(3))/2,dim(3)));
                        Xq(Xq<1)=1; Xq(Xq>dim(2))=dim(2);
                        Yq(Yq<1)=1; Yq(Yq>dim(1))=dim(1);
                        Zq(Zq<1)=1; Zq(Zq>dim(3))=dim(3);
                    end
                    if numel(shifts_temp) > 3
                        for dm = 1:3; shifts_up(:,:,:,dm) = interp3(shifts_temp(:,:,:,dm),Xq,Yq,Zq,'makima'); end
                        %tform = affine3d(diag([do([2,1,3])';1]));
                        %for dm = 1:3; shifts_up(:,:,:,dm) = imwarp(shifts_temp(:,:,:,dm),tform,'OutputView',imref3d([options.d1,options.d2,options.d3]),'SmoothEdges',true); end
                    else
                        for dm = 1:3; shifts_up(:,:,:,dm) = shifts_temp(dm); end
                    end
                    shifts_up(2:2:end,:,:,2) = shifts_up(2:2:end,:,:,2) + col_shift;
                    Mf = imwarp(Yt,-cat(4,shifts_up(:,:,:,2),shifts_up(:,:,:,1),shifts_up(:,:,:,3)),options.shifts_method,'FillValues',fill_value); 
                else
                    shifts_up = imresize(shifts_temp,[options.d1,options.d2]);
                    shifts_up(2:2:end,:,2) = shifts_up(2:2:end,:,2) + col_shift;
                    Mf = imwarp(Yt,-cat(3,shifts_up(:,:,2),shifts_up(:,:,1)),options.shifts_method,'FillValues',fill_value);  
                end    
             
             Mf(Mf<minY) = minY;
             Mf(Mf>maxY) = maxY;  
        end        
            
        if ~strcmpi(options.output_type,'mat')
            rem_mem = rem(t,options.mem_batch_size);
            if rem_mem == 0; rem_mem = options.mem_batch_size; end            
            if nd == 2; mem_buffer(:,:,rem_mem) = cast(Mf,data_type); end
            if nd == 3; mem_buffer(:,:,:,rem_mem) = cast(Mf,data_type); end
        end
        switch lower(options.output_type)
            case 'mat'
                if nd == 2; M_final(:,:,t) = cast(Mf,data_type); end
                if nd == 3; M_final(:,:,:,t) = cast(Mf,data_type); end
            case 'memmap'
                if rem_mem == options.mem_batch_size || t == T
                    if nd == 2; M_final.Y(:,:,t-rem_mem+1:t) = mem_buffer(:,:,1:rem_mem); end
                    if nd == 3; M_final.Y(:,:,:,t-rem_mem+1:t) = mem_buffer(:,:,:,1:rem_mem); end
                    M_final.Yr(:,t-rem_mem+1:t) = reshape(mem_buffer(1:d1*d2*d3*rem_mem),d1*d2*d3,rem_mem);
                end      
            case {'hdf5','h5'}
                if rem_mem == options.mem_batch_size || t == T
                    if nd == 2; h5write(options.h5_filename,['/',options.h5_groupname],mem_buffer(:,:,1:rem_mem),[ones(1,nd),t-rem_mem+1],[sizY(1:nd),rem_mem]); end
                    if nd == 3; h5write(options.h5_filename,['/',options.h5_groupname],mem_buffer(:,:,:,1:rem_mem),[ones(1,nd),t-rem_mem+1],[sizY(1:nd),rem_mem]); end
                end
            case {'tif','tiff'}
                if rem_mem == options.mem_batch_size || t == T
                    saveastiff(cast(mem_buffer(:,:,1:rem_mem),data_type),options.tiff_filename,opts_tiff);
                end                
        end         
        
        if mod(t,bin_width) == 0 && upd_template
            if print_msg
                str=[num2str(t), ' out of ', num2str(T), ' frames registered, iteration ', num2str(it), ' out of ', num2str(iter), '..'];
                refreshdisp(str, prevstr, t);
                prevstr=str; 
                %fprintf('%i out of %i frames registered, iteration %i out of %i \n',t,T,it,iter)
            end
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

if print_msg; fprintf('\n'); end
    
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
if print_msg; fprintf('done. \n'); end
end