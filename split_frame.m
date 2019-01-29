function Yf = split_frame(Y,patches)

df = diff(patches(1,:));
ndim = size(patches,2)/2;
Yf = zeros([df(1:2:end)+1,size(patches,1)],'like',Y);
otherdims = repmat({':'},1,ndims(Yf)-1+(size(patches,1)==1));

for i = 1:size(patches,1)
    indeces = cell(1,ndim);
    for j = 1:ndim
        indeces{j} = patches(i,2*j-1):patches(i,2*j);
    end
    Yf(otherdims{:},i) = Y(indeces{:});
end