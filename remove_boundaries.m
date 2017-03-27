function X = remove_boundaries(I,shifts,method,template,add_value)

if nargin < 5 || isempty(add_value);
    add_value = min(I(:));
end

if nargin < 3 || isempty(method);
    method = 'nan';
end

shifts = round(shifts);
sz = size(I);
nd = ndims(I);
X = I;

if strcmpi(method,'nan');
    method = 'zero';
    add_value = NaN;
end

switch lower(method)
    case 'zero'
        if shifts(1); X((1:abs(shifts(1)))*sign(shifts(1)) + (sz(1)+1)*(shifts(1)<0),:,:) = add_value; end
        if shifts(2); X(:,(1:abs(shifts(2)))*sign(shifts(2)) + (sz(2)+1)*(shifts(2)<0),:) = add_value; end
        if nd == 3 && shifts(3) 
            X(:,:,(1:abs(shifts(3)))*sign(shifts(3)) + (sz(3)+1)*(shifts(3)<0)) = add_value; 
        end        
    case 'copy'
        if shifts(1); X((1:abs(shifts(1)))*sign(shifts(1)) + (sz(1)+1)*(shifts(1)<0),:,:) = repmat(X((abs(shifts(1))+1)*sign(shifts(1)) + (sz(1)+1)*(shifts(1)<0),:,:),abs(shifts(1)),1,1); end
        if shifts(2); X(:,(1:abs(shifts(2)))*sign(shifts(2)) + (sz(2)+1)*(shifts(2)<0),:) = repmat(X(:,(abs(shifts(2))+1)*sign(shifts(2)) + (sz(2)+1)*(shifts(2)<0),:),1,abs(shifts(2)),1); end
        if nd == 3 && shifts(3) 
            X(:,:,(1:abs(shifts(3)))*sign(shifts(3)) + (sz(3)+1)*(shifts(3)<0)) = repmat(X(:,:,(abs(shifts(3))+1)*sign(shifts(3)) + (sz(3)+1)*(shifts(3)<0),:),1,1,abs(shifts(3))); 
        end
    case 'template'
        if shifts(1); X((1:abs(shifts(1)))*sign(shifts(1)) + (sz(1)+1)*(shifts(1)<0),:,:) = template((1:abs(shifts(1)))*sign(shifts(1)) + (sz(1)+1)*(shifts(1)<0),:,:); end
        if shifts(2); X(:,(1:abs(shifts(2)))*sign(shifts(2)) + (sz(2)+1)*(shifts(2)<0),:) = template(:,(1:abs(shifts(2)))*sign(shifts(2)) + (sz(2)+1)*(shifts(2)<0),:); end
        if nd == 3 && shifts(3) 
            X(:,:,(1:abs(shifts(3)))*sign(shifts(3)) + (sz(3)+1)*(shifts(3)<0)) = template(:,:,(1:abs(shifts(3)))*sign(shifts(3)) + (sz(3)+1)*(shifts(3)<0)); 
        end        
end