function [Ycon,ln] = concatenate_files(files,filename,output_type)

% concatenate a list of files into a single file
% INPUTS
% files:        list of a files as a struct array. Name of each file is in files(i).name
% filename:     output filename
% type:         output type ('mat' array loaded in RAM (default), 'hdf5' and 'tif' are currently supported)

% OUTPUTS
% Ycon:         concantenated file as an array if output_type == 'mat'
% ln:           length of each file (in frames)

if ~exist('output_type','var')
    output_type = 'mat';
end

if ~exist('filename','var') || isempty(filename)
    filename = ['concatenated_file.',output_type];
end

numFiles = length(files);
ln = zeros(numFiles,1);

Y1 = read_file(fullfile(files(1).folder, files(1).name));
sizY = size(Y1);
T1 = sizY(end);
data_type = class(Y1);
nd = max(ndims(Y1)-1,2);
sizY = sizY(1:nd);

switch lower(output_type)
    case 'mat'
        Ycon = cast([],data_type); %zeros([sizY,T1],data_type);       
    case {'hdf5','h5'}
        Ycon = filename;
        h5create(filename,'/mov',[sizY,Inf],'Chunksize',[sizY,1],'Datatype',data_type);
    case {'tif','tiff'}
        Ycon = filename;
        opts_tiff.append = true;
        opts_tiff.big = true;
    otherwise
        error('This filetype is currently not supported')
end   

for i = 1:numFiles
    Y_temp = read_file(fullfile(files(i).folder, files(i).name)); %cat(nd+1,Ycon,Y_temp);
    ln(i) = size(Y_temp,ndims(Y_temp));
    
    switch lower(output_type)
        case 'mat'
            Ycon = cat(nd+1,Ycon,Y_temp); 
        case {'hdf5','h5'}
            h5write(filename,'/mov',Y_temp,[ones(1,nd),sum(ln(1:i-1))+1],size(Y_temp));
        case {'tif','tiff'}
            saveastiff(cast(Y_temp,data_type),filename,opts_tiff);
    end         
end