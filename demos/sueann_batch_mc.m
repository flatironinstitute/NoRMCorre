clear;

foldername = '/Users/epnevmatikakis/Documents/Ca_datasets/Sueann/Sue/20161122';
files = dir([foldername,'/*.tif']);

%% set parameters
% options.grid_size = [128,128];             % size of patch in each direction
% options.bin_width = 50;                     % number of bins after which you update template
% options.mot_uf = 4;                         % upsampling factor for smaller patches
% options.us_fac = 20;                        % upsampling factor for subpixel registration
% options.method = {'median','mean'};         % averaging method for computing and updating templates
% options.overlap_pre = 32;                   % amount of overlap for each patch
% options.overlap_post = 16;                  % amount of overlap for each patch
% options.plot_flag = false;                  % flag for plotting results while correcting
% options.memmap = false;                     % save output in a .mat file
% options.iter = 1;
% options.max_shift = 15;
options = NoRMCorreSetParms('d1',512,'d2',512,'grid_size',128,'bin_width',50,'us_fac',20,'max_shift',15);
%% 
templates = cell(length(files),1);
%%
for i = 7:length(files) 
    if i > 7
        Y = (read_file([foldername,'/',files(i).name]));
    end
    cY = class(Y);
    if i == 6
        template = median(Y(:,:,1:50),3);
    end
    [M,shifts,template] = normcorre(Y,options,single(template));
    output_name = [foldername,'/',files(i).name(1:end-4),'_mc.tif'];
    %save_tiff(M,output_name,'single');
    saveastiff(int16(M), output_name);
    templates{i} = template;
    disp(i)
end

%%
tic;
data = memmap_file_sequence([foldername,'/mc']);
toc