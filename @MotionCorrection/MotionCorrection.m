classdef MotionCorrection < handle
    
    % MotionCorrection object for the NoRMCorre algorithm
    properties
        file ='';            % filename
        dims;                % dimensionality of the FOV
        nd;                  % number of dimensions (2d or 3d)
        T;                   % number of timesteps
        options;             % motion correction options structure
        shifts;              % inferred shifts
        shifts_x;            % shifts in x direction (matrix format)
        shifts_y;            % shifts in y direction (matrix format)
        shifts_z;            % shifts in z direction (matrix format)
        template;            % inferred template
        M;                   % motion corrected file
        meanY;               % mean image raw data
        meanM;               % mean image corrected data
        corrY                % correlation with mean raw data
        corrM;               % correlation with mean corrected data
        crispY;              % crispness of mean image (raw)
        crispM;              % crispness of mean image (corrected)
        col_shift = 0;       % relative shift due to bidirectional scanning
        gSig = 7;            % std of gaussian kernel for high pass filtering
        gSiz = 21;           % size of high pass filter
        file_orig;           % original file (used for 1p correction)        
    end
    
    methods                
        %% construct object, lazily load file and set dimensionality variables
        function obj = MotionCorrection(filename)
            obj.file = filename;
            if ischar(filename)
                Ytemp = read_file(filename,1,1);
                obj.dims = size(Ytemp);
            else    % filename is an array already loaded in memory
                sizY = size(filename);
                obj.T = sizY(end);          % assuming more than 1 frames
                obj.dims = sizY(1:end-1);                      
            end
            obj.nd = length(obj.dims);  
            obj.options.d1 = obj.dims(1);
            obj.options.d2 = obj.dims(2);
            if obj.nd == 3; obj.options.d3 = obj.dims(3); end
        end
        
        %% load and update options
        function optionsSet(obj, varargin)
            obj.options = NoRMCorreSetParms(obj.options,varargin{:});
        end
        
        %% motion correct serial
        function motionCorrectSerial(obj,options)
            if exist('options','var')
                obj.optionsSet(obj.options,options);
            end
            obj.options.d1 = obj.dims(1);
            obj.options.d2 = obj.dims(2);
            if obj.nd == 3; obj.options.d3 = obj.dims(3); end
            [obj.M,obj.shifts,obj.template,obj.options,obj.col_shift] = normcorre(obj.file,obj.options);
            obj.T = length(obj.shifts);
            shifts_res = cat(ndims(obj.shifts(1).shifts)+1,obj.shifts(:).shifts);
            shifts_res = reshape(shifts_res,[],obj.nd,obj.T);
            obj.shifts_x = squeeze(shifts_res(:,1,:))';
            obj.shifts_y = squeeze(shifts_res(:,2,:))';
            if obj.nd == 3; obj.shifts_z = squeeze(shifts_res(:,3,:))'; end
        end
        
        %% motion correct Parallel(obj)
        function motionCorrectParallel(obj,options)
            if exist('options','var')
                obj.optionsSet(obj.options,options);
            end
            obj.options.d1 = obj.dims(1);
            obj.options.d2 = obj.dims(2);
            if obj.nd == 3; obj.options.d3 = obj.dims(3); end
            [obj.M,obj.shifts,obj.template,obj.options,obj.col_shift] = normcorre_batch(obj.file,obj.options);
            obj.T = length(obj.shifts);
            shifts_res = cat(ndims(obj.shifts(1).shifts)+1,obj.shifts(:).shifts);
            shifts_res = reshape(shifts_res,[],obj.nd,obj.T);
            obj.shifts_x = squeeze(shifts_res(:,1,:))';
            obj.shifts_y = squeeze(shifts_res(:,2,:))';
            if obj.nd == 3; obj.shifts_z = squeeze(shifts_res(:,3,:))'; end
        end
        
        %% apply shifts to a different file
        function M = applyShifts(obj,newfile)
            M = apply_shifts(newfile,obj.shifts,obj.options,0,0,0,obj.col_shift);
        end
        
        %% high pass filter
        function HPF(obj,gSig,gSiz)
            if exist('gSig','var'); obj.gSig = gSig; end
            if exist('gSiz','var'); obj.gSiz = gSiz; end
            psf = fspecial('gaussian', round(gSiz), gSig);
            ind_nonzero = (psf(:)>=max(psf(:,1)));
            psf = psf-mean(psf(ind_nonzero));
            psf(~ind_nonzero) = 0;               
            if ischar(obj.file)
                obj.file_orig = read_file(obj.file);
            else
                obj.file_orig = obj.file;                
            end
            obj.file = imfilter(single(obj.file_orig),psf,'symmetric');
        end
        
        %% compute mean
        function computeMean(obj)
            if obj.nd == 2; batch = 1000; elseif obj.nd == 3; batch = 10; end
            if isempty(obj.meanY)
                if ischar(obj.file)
                    obj.meanY = zeros(obj.dims,'single');
                    for t = 1:batch:obj.T
                        Y_temp = single(read_file(obj.file,t,batch));
                        obj.meanY = ((t-1)*obj.meanY+size(Y_temp,obj.nd+1)*mean(Y_temp,obj.nd+1))...
                                        /(t-1+size(Y_temp,obj.nd+1));
                    end
                else
                    obj.meanY = mean(obj.Y,obj.nd+1);
                end
            end
            if isempty(obj.meanM)
                if ischar(obj.M)
                    obj.meanM = zeros(obj.dims,'single');
                    for t = 1:batch:obj.T
                        M_temp = single(read_file(obj.M,t,batch));
                        obj.meanM = ((t-1)*obj.meanM+size(M_temp,obj.nd+1)*mean(M_temp,obj.nd+1))...
                                        /(t-1+size(M_temp,obj.nd+1));
                    end
                else
                    obj.meanM = mean(obj.M,obj.nd+1);
                end
            end
        end
        
        %% crispness (norm of gradient of mean)
        function crispness(obj,bnd)
            if ~exist('bnd','var'); bnd = zeros(6,1); end
            if isscalar(bnd); bnd = bnd*ones(6,1); end
            computeMean(obj);
            if obj.nd == 2
                mean_Y = obj.meanY(bnd(1)+1:end-bnd(2),bnd(3)+1:end-bnd(4));
                mean_M = obj.meanM(bnd(1)+1:end-bnd(2),bnd(3)+1:end-bnd(4));
            elseif obj.nd == 3
                mean_Y = obj.meanY(bnd(1)+1:end-bnd(2),bnd(3)+1:end-bnd(4),bnd(5)+1:end-bnd(6));
                mean_M = obj.meanM(bnd(1)+1:end-bnd(2),bnd(3)+1:end-bnd(4),bnd(5)+1:end-bnd(6));
            end
            if ismatrix(obj.meanY)
                [gx,gy] = gradient(mean_Y); 
                ng = norm(sqrt(gx.^2+gy.^2),'fro');
            elseif obj.nd == 3        
                [gx,gy,gz] = gradient(mean_Y); 
                ng = sqrt(sum(gx(:).^2+gy(:).^2+gz(:).^2));
            end
            obj.crispY = ng;
            if ismatrix(obj.meanM)
                [gx,gy] = gradient(mean_M); 
                ng = norm(sqrt(gx.^2+gy.^2),'fro');
            elseif obj.nd == 3
                [gx,gy,gz] = gradient(mean_M); 
                ng = sqrt(sum(gx(:).^2+gy(:).^2+gz(:).^2));
            end
            obj.crispM = ng;
        end
        
        %% correlation with mean
        function correlationMean(obj,bnd)
            if ~exist('bnd','var'); bnd = zeros(6,1); end
            if isscalar(bnd); bnd = bnd*ones(6,1); end
            computeMean(obj);
            if obj.nd == 2
                mean_Y = obj.meanY(bnd(1)+1:end-bnd(2),bnd(3)+1:end-bnd(4));
                mean_M = obj.meanM(bnd(1)+1:end-bnd(2),bnd(3)+1:end-bnd(4));
            elseif obj.nd == 3
                mean_Y = obj.meanY(bnd(1)+1:end-bnd(2),bnd(3)+1:end-bnd(4),bnd(5)+1:end-bnd(6));
                mean_M = obj.meanM(bnd(1)+1:end-bnd(2),bnd(3)+1:end-bnd(4),bnd(5)+1:end-bnd(6));
            end
            if obj.nd == 2; batch = 1000; elseif obj.nd == 3; batch = 10; end
            if ~isempty(obj.file)
                if ischar(obj.file)
                    obj.corrY = zeros(obj.T,1,'single');
                    for t = 1:batch:obj.T
                        Y_temp = single(read_file(obj.file,t,batch));
                        if obj.nd == 2; Y_temp = Y_temp(bnd(1)+1:end-bnd(2),bnd(3)+1:end-bnd(4),:); end
                        if obj.nd == 3; Y_temp = Y_temp(bnd(1)+1:end-bnd(2),bnd(3)+1:end-bnd(4),bnd(5)+1:end-bnd(6),:); end
                        obj.corrY(t:t+size(Y_temp,obj.nd+1)-1) = corr(reshape(Y_temp,[],size(Y_temp,obj.nd+1)),mean_Y(:));
                    end
                else
                    if obj.nd == 2; obj.corrY = corr(mean_Y(:),single(reshape(obj.Y(bnd(1)+1:end-bnd(2),bnd(3)+1:end-bnd(4),:),[],obj.T))); end
                    if obj.nd == 3; obj.corrY = corr(mean_Y(:),single(reshape(obj.Y(bnd(1)+1:end-bnd(2),bnd(3)+1:end-bnd(4),bnd(5)+1:end-bnd(6),:),[],obj.T))); end
                end
            end
            if ~isempty(obj.M)
                if ischar(obj.M)
                    obj.meanM = zeros(obj.T,1,'single');
                    for t = 1:batch:obj.T
                        M_temp = single(read_file(obj.M,t,batch));
                        if obj.nd == 2; M_temp = M_temp(bnd(1)+1:end-bnd(2),bnd(3)+1:end-bnd(4),:); end
                        if obj.nd == 3; M_temp = M_temp(bnd(1)+1:end-bnd(2),bnd(3)+1:end-bnd(4),bnd(5)+1:end-bnd(6),:); end
                        obj.corrM(t:t+size(M_temp,obj.nd+1)-1) = corr(reshape(M_temp,[],size(M_temp,obj.nd+1)),mean_M(:));
                    end
                else
                    if obj.nd == 2; obj.corrM = corr(mean_M(:),single(reshape(obj.M(bnd(1)+1:end-bnd(2),bnd(3)+1:end-bnd(4),:),[],obj.T))); end
                    if obj.nd == 3; obj.corrM = corr(mean_M(:),single(reshape(obj.M(bnd(1)+1:end-bnd(2),bnd(3)+1:end-bnd(4),bnd(5)+1:end-bnd(6),:),[],obj.T))); end
                end
            end
        end
        
    end
end