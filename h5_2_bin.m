function h5_2_bin(path_to_file,output_name,precision,batch_size)

info = hdf5info(path_to_file);
dims = info.GroupHierarchy.Datasets.Dims;
T = dims(end);
if nargin < 4 || isempty(batch_size);
    batch_size = T;
else
    batch_size = min(batch_size,T);
end
if nargin < 3 || isempty(precision)
    precision = 'uint16';
end

fid = fopen(output_name,'w+');
fclose(fid);
for t = 1:batch_size:T
    fid = fopen(output_name,'a');
    num2read = min(batch_size,T-t+1);
    imData = h5read(path_to_file,'/mov',[ones(1,length(dims)-1),t],[dims(1:end-1),num2read]);
    imData = cast(imData,precision);
    fseek(fid,0,'eof');
    fwrite(fid,imData,precision);
    fclose(fid);
end
%fclose(fid);