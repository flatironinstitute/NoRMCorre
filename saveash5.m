function saveash5(data, path_to_file, options)

% save the variable stored in data in a h5 file named path_to_file (includes the
% path to the directory). The file is saved with the same class as data.
% options.groupname: the name of the group (default: /mov)
% options.append: append to existing file if it exists (default: True)

if ~exist('path_to_file','var'); path_to_file = 'data.h5'; end
[pathstr, name, ~] = fileparts(path_to_file);
path_to_file = fullfile(pathstr,[name,'.h5']); % make sure file has an h5 extension
data_type = class(data);

sizY = size(data);
nd = ndims(data)-1;

defoptions.groupname = '/mov';
defoptions.append = true;

if ~exist('options','var'); options = defoptions; end
if ~isfield(options,'groupname'); options.groupname = defoptions.groupname; end
if ~isfield(options,'groupname'); options.append = defoptions.append; end

if ~exist(path_to_file,'file') || ~options.append
    if ~options.append && exist('path_to_file','file')
        warning('Rewriting existing file')
    end
    h5create(path_to_file,options.groupname,[sizY(1:nd),Inf],'Chunksize',[sizY(1:nd),100],'Datatype',data_type);
    start_point = ones(1,nd+1);
else
    fprintf('Appending to existing file \n')
    info = h5info(path_to_file);
    options.groupname = ['/',info.Datasets.Name];
    start_point = [ones(1,nd),info.Datasets.Dataspace.Size(end)+1];
end

h5write(path_to_file,options.groupname,data,start_point,size(data));