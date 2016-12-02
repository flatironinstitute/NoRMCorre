name = '/Users/epnevmatikakis/Dropbox (Simons Foundation)/fly1_runD_malefly_10_5_10/skewed_1.tiff';

Z = bigread2(name);
Z2 = Z(:,511:end,:);

Y = zeros([size(Z2),250],'uint16');
Y(:,:,:,1) = Z2;

%%

for t = 2:300
    name = ['/Users/epnevmatikakis/Dropbox (Simons Foundation)/fly1_runD_malefly_10_5_10/skewed_',num2str(t),'.tiff'];
    Z = bigread2(name);
    Y(:,:,:,t) = Z(:,511:end,:);
    disp(t);
end

tic; savefast('evan2_trun_full.mat','Y'); toc
data = matfile('evan2_trun.mat','Writable',true);
%%

options.grid_size = [64,64,56];             % size of patch in each direction
options.bin_width = 50;                     % number of bins after which you update template
options.mot_uf = 2;                         % upsampling factor for smaller patches
options.us_fac = 4;                         % upsampling factor for subpixel registration
options.method = {'median','mean'};         % averaging method for computing and updating templates
options.overlap_pre = 12;                   % amount of overlap for each patch
options.iter = 1;                           % number of passes (set to 1 for one-pass online processing)
options.plot_flag = false;                   % flag for plotting results while correcting
options.memmap = true;                     % save output in a .mat file
options.make_avi = false;                   % make a movie showing the results
options.fr = 15;                            % frame rate for movie
options.filename = 'evan2_nr2_new.mat';
options.use_parallel = false;
%options.name = [name,'_corrected.avi'];     % name for movie

%% perform motion correction
tt1 = tic;
%profile on;
[M_final,shifts,template] = online_motion_correction_patches(data,options);
%profile off;
%profile viewer
toc(tt1)

%%
[cY,mY,vY] = motion_metrics(data,[10,10,0]);
%%
[cM,mM,vM] = motion_metrics(M_final,[10,10,0]);