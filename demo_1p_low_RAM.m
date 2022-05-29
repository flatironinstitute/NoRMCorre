% demo file for applying the NoRMCorre motion correction algorithm on 
% 1-photon widefield imaging data using low memory (good for long datasets)
% Example file is provided from the miniscope project page
% www.miniscope.org

clear;
gcp;
%% read data and convert to double
name = 'msCam13.avi';
if ~exist(name,'file')  % download file if it doesn't exist in the directory
    url = 'https://caiman.flatironinstitute.org/~neuro/normcorre_datasets/msCam13.avi';
    fprintf('downloading the file...');
    outfilename = websave(name,url);
    fprintf('done.');
end
frame = read_file(name,1,1);
[d1,d2] = size(frame);

%% perform some sort of deblurring/high pass filtering
% The function does not load the whole file in memory. Instead it loads 
% chunks of the file and then saves the high pass filtered version in a 
% h5 file.

gSig = 7; 
gSiz = 3*gSig; 
psf = fspecial('gaussian', round(2*gSiz), gSig);
ind_nonzero = (psf(:)>=max(psf(:,1)));
psf = psf-mean(psf(ind_nonzero));
psf(~ind_nonzero) = 0;   % only use pixels within the center disk

[filepath,file_name,ext] = fileparts(name);
h5_name = fullfile(filepath,[file_name,'_filtered_data.h5']);
chunksize = 750;    % read 500 frames at a time
cnt = 1;
while (1)  % read filter and save file in chunks
    Yf = single(read_file(name,cnt,chunksize));
    if isempty(Yf)
        break
    else
        Y = imfilter(Yf,psf,'symmetric');
        saveash5(Y,h5_name);
        cnt = cnt + size(Y,ndims(Y));
    end
    disp(cnt)
end
%% first try out rigid motion correction
    % exclude boundaries due to high pass filtering effects
options_r = NoRMCorreSetParms('d1',d1,'d2',d2,'bin_width',200,'max_shift',20,'iter',1,'correct_bidir',false);

%% register using the high pass filtered data and apply shifts to original data
tic; [M1,shifts1,template1] = normcorre_batch(h5_name,options_r); toc % register filtered data
    % exclude boundaries due to high pass filtering effects
    
% if you save the file directly in memory make sure you save it with a 
% name that does not exist. Change options_r.tiff_filename 
% or options_r.h5_filename accordingly.

tic; Mr = apply_shifts(name,shifts1,options_r); toc % apply shifts to full dataset
    
% you can only save the motion corrected file directly in memory by
% setting options_r.output_type = 'tiff' or 'h5' and selecting an
% appropriate name through options_r.tiff_filename or options_r.h5_filename