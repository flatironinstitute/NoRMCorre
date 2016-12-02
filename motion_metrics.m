function [cY,mY,vY] = motion_metrics(Y,bnd,var_name)

if nargin == 1 || isempty(bnd); bnd = zeros(6,1); end
if nargin == 2 || isempty(batch_size); batch_size = 1; end

memmap = isobject(Y);

if memmap
    if ~exist('var_name','var')
        try sizY = size(Y,'Z'); var_name = 'Y'; catch; sizY = size(Y,'M'); var_name = 'M'; end
    end
else
    sizY = size(Y);
end

dimsY = length(sizY);
d = prod(sizY(1:end-1));
T = sizY(end);

if isscalar(bnd); bnd = ones(2*(dimsY-1),1)*bnd; end
if dimsY == 3; sizY(3) = 1; bnd(5:6) = 0; end

if memmap
    cY = zeros(T,1);
    mY = zeros(sizY(1:end-1));

    if nargout == 3
        vY = zeros(sizY(1:end-1));
    end

    for t = 1:T
        y_temp = single(load_data(Y,t));
        delta = y_temp - mY;
        mY = mY + delta/t;
        if nargout == 3
            vY = vY + delta.*(y_temp - mY);
        end    
    end

    for t = 1:T
        y_temp = single(load_data(Y,t));
        y_temp = y_temp(bnd(1)+1:sizY(1)-bnd(2),bnd(3)+1:sizY(2)-bnd(4),bnd(5)+1:sizY(3)-bnd(6));
        m_temp = mY(bnd(1)+1:sizY(1)-bnd(2),bnd(3)+1:sizY(2)-bnd(4),bnd(5)+1:sizY(3)-bnd(6));
        cY(t) = corr(y_temp(:),m_temp(:));
    end
else
    Y = single(Y);
    nd = ndims(Y)-1;
    mY = mean(Y,nd+1);
    if nd == 2
        Yr = Y(bnd(1)+1:sizY(1)-bnd(2),bnd(3)+1:sizY(2)-bnd(4),:);
        mYr = mean(Yr,nd+1); %mY(bnd(1)+1:sizY(1)-bnd(2),bnd(3)+1:sizY(2)-bnd(4));
    else
        Yr = Y(bnd(1)+1:sizY(1)-bnd(2),bnd(3)+1:sizY(2)-bnd(4),bnd(5)+1:sizY(3)-bnd(6),:);
        mYr = mY(bnd(1)+1:sizY(1)-bnd(2),bnd(3)+1:sizY(2)-bnd(4),bnd(5)+1:sizY(3)-bnd(6));
    end
    cY = corr(reshape(Yr,[],T),reshape(mYr,[],1));
    if nargout == 3; vY = var(Y,[],nd+1)*(T-1); end
end
    
if nargout == 3; vY = vY/T; end

    function y_temp = load_data(X,t)
        if ~memmap
            y_temp = reshape(X(d*(t-1) + (1:d)),sizY(1:end-1));
        else
            if strcmp(var_name,'Y')
                if dimsY == 3
                    y_temp = double(X.Y(:,:,t));
                else
                    y_temp = double(X.Y(:,:,:,t));
                end
            elseif strcmp(var_name,'M')
                if dimsY == 3
                    y_temp = double(X.M(:,:,t));
                else
                    y_temp = double(X.M(:,:,:,t));
                end
            else
                error('unknown variable name')
            end
        end
    end
end