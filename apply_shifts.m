function M_final = apply_shifts(Y,shifts,options,td1,td2,td3,col_shift)

% apply shifts using fft/cubic or linear interpolation

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

if isa(Y,'char')
    [~,~,ext] = fileparts(Y);
    ext = ext(2:end);
    if strcmpi(ext,'tif') || strcmpi(ext,'tiff')
        tiffInfo = imfinfo(Y);
        sizY = [tiffInfo(1).Height,tiffInfo(1).Width,length(tiffInfo)];
        filetype = 'tif';
        data_type = class(imread(Y,'Index',1,'Info',tiffInfo));
    elseif strcmpi(ext,'mat')
        filetype = 'mem';
        Y = matfile(Y,'Writable',true);
        sizY = size(Y,'Y');
        details = whos(Y,'Y');
        data_type = details.class;
    elseif strcmpi(ext,'hdf5') || strcmpi(ext,'h5')
        filetype = 'hdf5';
        fileinfo = hdf5info(Y);
        data_name = fileinfo.GroupHierarchy.Datasets.Name;
        sizY = fileinfo.GroupHierarchy.Datasets.Dims;
        data_type = class(read_file(Y,1,1));
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
        data_type = 'single';
    elseif strcmpi(ext,'avi')
        filetype = 'avi';
        v = VideoReader(Y);
        T = v.Duration*v.FrameRate;
        sizY = [v.Height,v.Width,T];
        data_type = class(readFrame(v));
    end    
elseif isobject(Y)
    filetype = 'mem';
    sizY = size(Y,'Y');
    details = whos(Y,'Y');
    data_type = details.class;
else % array loaded in memory
    filetype = 'mat';
    data_type = class(Y);
    Y = single(Y);
    sizY = size(Y);
end
if strcmpi(options.boundary,'nan')
    fill_value = NaN;
else
    fill_value = options.add_value;
end

T = length(shifts);

if sizY(end) == T && T > 1
    flag_constant = false;
    nd = length(sizY)-1;   
    sizY = sizY(1:end-1);
else
    flag_constant = true;
    nd = length(sizY);
end
d1 = sizY(1); d2 = sizY(2);
if nd == 2; d3 = 1; else d3 = sizY(3); end

if ~isfield(options, 'print_msg') || isempty(options.print_msg)
    print_msg = true;
else
    print_msg = options.print_msg;
end

if strcmpi(options.shifts_method,'fft')
    % precompute some quantities that are used repetitively for template matching and applying shifts

    [xx_s,xx_f,yy_s,yy_f,zz_s,zz_f,xx_us,xx_uf,yy_us,yy_uf,zz_us,zz_uf] = construct_grid(options.grid_size,options.mot_uf,options.d1,options.d2,options.d3,options.min_patch_size);
    xx_us = xx_us + td1; xx_us(1) = 1;
    yy_us = yy_us + td2; yy_us(1) = 1;
    zz_us = zz_us + td3; zz_us(1) = 1;
    xx_uf = xx_uf + td1; xx_uf(end) = d1;
    yy_uf = yy_uf + td2; yy_uf(end) = d2;
    zz_uf = zz_uf + td3; zz_uf(end) = d3;

    temp_cell = mat2cell_ov(zeros(d1,d2,d3,'single'),xx_us,xx_uf,yy_us,yy_uf,zz_us,zz_uf,options.overlap_post,[d1,d2,d3]);
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
    shift_fun = @(yfft,shfts,ph,nr,nc,np) shift_reconstruct(yfft,shfts,ph,options.us_fac,nr,nc,np,options.boundary,0);
    Xq = []; Yq = []; Zq = [];
else
    if nd == 3
        dim = [d1,d2,d3];
        ds = size(shifts(1).shifts);
        do = [d1,d2,d3,1]./size(shifts(1).shifts);
        %tform = affine3d(diag([do([2,1,3])';1]));
        [Xq,Yq,Zq] = meshgrid(linspace((1+1/do(2))/2,ds(2)+(1-1/do(2))/2,dim(2)),linspace((1+1/do(1))/2,ds(1)+(1-1/do(1))/2,dim(1)),linspace((1+1/do(3))/2,ds(3)+(1-1/do(3))/2,dim(3)));
    else
        Xq = []; Yq = []; Zq = [];
    end
end

switch lower(options.output_type)
    case 'mat'
        M_final = zeros([sizY(1:nd),T],data_type);
    case 'memmap'
        M_final = matfile(options.mem_filename,'Writable',true);
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


if exist('col_shift','var'); options.col_shift = col_shift; end
if ~isempty(options.col_shift) 
    col_shift = options.col_shift; 
    options.correct_bidir = false; 
elseif ~options.correct_bidir
    col_shift = 0;
end
if options.correct_bidir
    col_shift = correct_bidirectional_offset(Y,options.nFrames,options.bidir_us);
end

if col_shift
    if strcmpi(options.shifts_method,'fft')
        options.shifts_method = 'cubic';
        if print_msg; fprintf('Offset %1.1f pixels due to bidirectional scanning detected. Cubic shifts will be applied. \n',col_shift); end
    end
end

if print_msg; prevstr = []; end
bin_width = min([options.mem_batch_size,T,ceil((512^2*3000)/(d1*d2*d3))]);
for t = 1:bin_width:T
    switch filetype
        case 'tif'
            Ytm = single(read_file(Y, t, min(t+bin_width-1,T)-t+1, [], tiffInfo));
%             Ytm = zeros(sizY(1),sizY(2),min(t+bin_width-1,T)-t+1,'single');
%             for tt = 1:min(t+bin_width-1,T)-t+1
%                 Ytm(:,:,tt) = single(imread(Y,'Index',t+tt-1,'Info',tiffInfo));
%             end
        case 'hdf5'
            Ytm = single(h5read(Y,data_name,[ones(1,length(sizY)-1),t],[sizY(1:end-1),min(t+bin_width-1,T)-t+1]));
        case 'mem'
            if nd == 2; Ytm = single(Y.Y(:,:,t:min(t+bin_width-1,T))); end
            if nd == 3; Ytm = single(Y.Y(:,:,:,t:min(t+bin_width-1,T))); end
        case 'mat'
            if nd == 2; Ytm = single(Y(:,:,t:min(t+bin_width-1,T))); end
            if nd == 3; Ytm = single(Y(:,:,:,t:min(t+bin_width-1,T))); end
        case 'raw'
            Ytm = read_raw_file(Y,t,min(t+bin_width-1,T)-t+1,FOV,bitsize);
        case 'avi'
            Ytm = read_file(Y,t,min(t+bin_width-1,T)-t+1);
    end
%    if ~flag_constant
    if nd == 2; Ytc = mat2cell(Ytm,d1,d2,ones(1,size(Ytm,3))); end
    if nd == 3; Ytc = mat2cell(Ytm,d1,d2,d3,ones(1,size(Ytm,4))); end
    
    Mf = cell(size(Ytc));
    lY = length(Ytc);
    shifts_temp = shifts(t:t+lY-1);
    
    switch lower(options.shifts_method)
        case 'fft'            
            parfor ii = 1:lY 
                Yc = mat2cell_ov(Ytc{ii},xx_us,xx_uf,yy_us,yy_uf,zz_us,zz_uf,options.overlap_post,[d1,d2,d3]);
                Yfft = cellfun(@(x) fftn(x),Yc,'un',0);
                minY = min(Ytc{ii}(:));
                maxY = max(Ytc{ii}(:));

                if all([length(xx_s),length(yy_s),length(zz_s)] == 1)
                    M_fin = shift_reconstruct(Yfft{1},shifts_temp(ii).shifts,shifts_temp(ii).diff,options.us_fac,Nr{1},Nc{1},Np{1},options.boundary,0);
                    Mf{ii} = M_fin;
                else
                    shifts_up = shifts_temp(ii).shifts_up;
                    shifts_cell = mat2cell(shifts_up,ones(length(xx_uf),1),ones(length(yy_uf),1),ones(length(zz_uf),1),nd);
                    diff_cell = num2cell(shifts_temp(ii).diff);
                    M_fin = cellfun(shift_fun,Yfft,shifts_cell,diff_cell,Nr,Nc,Np,'un',0);

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
                Mf{ii}(Mf{ii}<minY) = minY;
                Mf{ii}(Mf{ii}>maxY) = maxY;    
            end
        otherwise
            parfor ii = 1:lY
                minY = min(Ytc{ii}(:));
                maxY = max(Ytc{ii}(:));
                shifts_temp(ii).shifts_up = shifts_temp(ii).shifts;
                if nd == 3                                    
                    shifts_up = zeros([d1,d2,d3,3]);
                    if numel(shifts_temp(ii).shifts) > 3
                        %tform = affine3d(diag([options.mot_uf(:);1]));
                        %tform = affine3d(diag([do([2,1,3])';1]));
                        %for dm = 1:3; shifts_up(:,:,:,dm) = imwarp(shifts_temp(ii).shifts(:,:,:,dm),tform,'OutputView',imref3d([d1,d2,d3]),'SmoothEdges',true); end
                        for dm = 1:3; shifts_up(:,:,:,dm) = interp3(shifts_temp(ii).shifts(:,:,:,dm),Xq,Yq,Zq,'makima'); end
                    else
                        for dm = 1:3; shifts_up(:,:,:,dm) = shifts_temp(ii).shifts(:,:,:,dm); end
                    end
                    shifts_up(2:2:end,:,:,2) = shifts_up(2:2:end,:,:,2) + col_shift;
                    Mf{ii} = imwarp(Ytc{ii},-cat(4,shifts_up(:,:,:,2),shifts_up(:,:,:,1),shifts_up(:,:,:,3)),options.shifts_method,'FillValues',fill_value);
                else
                    shifts_up = imresize(shifts_temp(ii).shifts,[options.d1,options.d2]);
                    shifts_up(2:2:end,:,2) = shifts_up(2:2:end,:,2) + col_shift;
                    Mf{ii} = imwarp(Ytc{ii},-cat(3,shifts_up(:,:,2),shifts_up(:,:,1)),options.shifts_method,'FillValues',fill_value);  
                end
                Mf{ii}(Mf{ii}<minY) = minY;
                Mf{ii}(Mf{ii}>maxY) = maxY;
            end
    end
    
    Mf = cast(cell2mat(Mf),data_type);
    
    switch lower(options.output_type)
        case 'mat'
            if nd == 2; M_final(:,:,t:min(t+bin_width-1,T)) = Mf; end
            if nd == 3; M_final(:,:,:,t:min(t+bin_width-1,T)) = Mf; end
        case 'memmap'
            if nd == 2; M_final.Y(:,:,t:min(t+bin_width-1,T)) = Mf; end
            if nd == 3; M_final.Y(:,:,:,t:min(t+bin_width-1,T)) = Mf; end
            M_final.Yr(:,t:min(t+bin_width-1,T)) = reshape(Mf,d1*d2*d3,[]);
        case {'hdf5','h5'}
            rem_mem = min(bin_width,T-t+1);
            if nd == 2; h5write(options.h5_filename,['/',options.h5_groupname],Mf,[ones(1,nd),t],[sizY(1:nd),rem_mem]); end
            if nd == 3; h5write(options.h5_filename,['/',options.h5_groupname],Mf,[ones(1,nd),t],[sizY(1:nd),rem_mem]); end
        case {'tif','tiff'}
            saveastiff(cast(Mf,data_type),options.tiff_filename,opts_tiff);
    end
    
    if print_msg
        str = sprintf('%i out of %i frames registered \n',t+lY-1,T);
        refreshdisp(str, prevstr, t);
        prevstr=str;
    end
end

if print_msg; fprintf('\n'); end
