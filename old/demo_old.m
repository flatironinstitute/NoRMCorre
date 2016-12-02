clear
%name = '/Users/epnevmatikakis/Documents/Python/github/Constrained_NMF/PPC.tif';
name = '/Users/epnevmatikakis/Downloads/efficient_subpixel_registration/J115_2015-12-09_L01_074.tif'; % set path to file
%name = '/Users/epnevmatikakis/downloads/efficient_subpixel_registration/M_FLUO_t.tif'; % set path to file
%name = 'M_FLUO_4.tif';
%name = '/Users/epnevmatikakis/Documents/Ca_datasets/Darcy/Single_150um_024.tif';
Y = bigread2(name); % read the file (optional, you can also pass the path in the function instead of Y)
Y = double(Y);      % convert to double precision 

%% set parameters

options.grid_size = [128,128];              % size of patch in each direction
options.bin_width = 20;                     % number of bins after which you update template
options.mot_uf = 4;                         % upsampling factor for smaller patches
options.us_fac = 6;                         % upsampling factor for subpixel registration
options.method = {'mean','mean'};           % averaging method for computing and updating templates
options.overlap_pre = 8;                    % amount of overlap for each patch
options.iter = 1;                           % number of passes (set to 1 for one-pass online processing)
options.plot_flag = true;                  % flag for plotting results while correcting
options.memmap = false;                     % save output in a .mat file
options.make_avi = false;                   % make a movie showing the results
options.fr = 15;                            % frame rate for movie
options.use_parallel = false;
options.name = [name,'_corrected.avi'];     % name for movie
options.iter = 1;
options.max_shift = 10;
%% perform motion correction
%tic;
profile on;
%[M_final4,shifts4,template4,shifts_up4] = online_motion_correction_patches(Y,options);
tic; [M_final1,shifts1,template1] = online_motion_correction_patches_global(Y,options); toc
%tic; [M_final2,shifts2,template2,shifts_up2] = online_motion_correction_patches(Y,options); toc

profile off;
profile viewer
%toc
%% plot data

nY = quantile(Y(:),0.005);
mY = quantile(Y(:),0.995);

figure;
for t = 1:2:T
    subplot(121);imagesc(Y(:,:,t),[nY,mY]);
    title(sprintf('Frame %i out of %i',t,T),'fontweight','bold','fontsize',14); colormap('bone')
    subplot(122);imagesc(M_final(:,:,t),[nY,mY]);
    drawnow;
    pause(0.05);
end