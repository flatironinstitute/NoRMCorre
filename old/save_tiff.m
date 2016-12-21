function save_tiff(Y,output_name,type)

type = lower(type);

switch type
    case 'uint16'
        Y = uint16(Y);
    case 'uint32'
        Y = uint32(Y);
    case 'int16'
        Y = int16(Y);
    case 'int32'
        Y = int32(Y);
    case 'single'
        Y = single(Y);
    otherwise
        Y = double(Y);
end

sizY = size(Y);
nd = ndims(Y);
T = sizY(nd);
d = prod(sizY(1:nd-1));
for t = 1:T   
    imwrite(reshape(Y((t-1)*d+(1:d)),sizY(1:nd-1)), output_name, 'WriteMode', 'append','Compression','none');
end