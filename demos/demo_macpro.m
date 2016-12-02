clear
%name = 'M_FLUO_4.tif';
%name = '/Users/epnevmatikakis/Documents/Ca_datasets/SIMA/2D_example.tif';
%name = '/Users/epnevmatikakis/Documents/Python/github/Constrained_NMF/PPC.tif';
name = '/Users/epnevmatikakis/Documents/Ca_datasets/Sueann/k56_20160608_RSM_125um_41mW_zoom2p2_00001_00034.tif';
Y = bigread2(name); % read the file (optional, you can also pass the path in the function instead of Y)
Y = double(Y);      % convert to double precision 
T = size(Y,ndims(Y));
%% set parameters (first try out rigid motion correction)

options.grid_size = [size(Y,1),size(Y,2)];  % size of patch in each direction
options.bin_width = 30;                     % number of bins after which you update template
options.mot_uf = 2;                         % upsampling factor for smaller patches
options.us_fac = 10;                         % upsampling factor for subpixel registration
options.method = {'median','mean'};         % averaging method for computing and updating templates
options.overlap_pre = 16;                   % amount of overlap for each patch
options.overlap_post = 16;                  % amount of overlap for each patch
options.iter = 2;                           % number of passes (set to 1 for one-pass online processing)
options.plot_flag = false;                  % flag for plotting results while correcting
options.memmap = false;                     % save output in a .mat file
options.iter = 1;
options.max_shift = 15;

%% perform motion correction
tic; [M1,shifts1,template1] = online_motion_correction_patches(Y,options); toc

%% now try non-rigid
options.grid_size = [32,128];                % size of patch in each direction
options.mot_uf = 4;                         % further upsample by a given factor
tic; [M2,shifts3,template3] = online_motion_correction_patches(Y,options); toc
options.grid_size = [16,16];                % size of patch in each direction
options.plot_flag = false;                  % flag for plotting results while correcting

options.mot_uf = 2;                         % further upsample by a given factor
profile on
tic; [M4,shifts4,template4] = online_motion_correction_patches(Y,options); toc
profile off
profile viewer
% These results are with only one pass of the data completely online. In
% practice we can improve if we do multiple passes (with options.iter > 1)
% or do fancy things like, first correct with rigid which is more robust to
% large perturbation, then refine by correcting with non-rigid
%% plot data

nnY = quantile(Y(:),0.005);
mmY = quantile(Y(:),0.995);

[cY,mY,vY] = motion_metrics(Y,10);
[cM1,mM1,vM1] = motion_metrics(M1,10);
[cM2,mM2,vM2] = motion_metrics(M2,10);
[cM3,mM3,vM3] = motion_metrics(M3,10);

T = length(cY);
%% plot metrics
figure;
    subplot(2,3,1); imagesc(mY,[nnY,mmY]);  axis equal; axis tight; axis off; title('mean raw data','fontsize',14,'fontweight','bold')
    subplot(2,3,2); imagesc(mM1,[nnY,mmY]);  axis equal; axis tight; axis off; title('mean rigid corrected','fontsize',14,'fontweight','bold')
    subplot(2,3,3); imagesc(mM2,[nnY,mmY]); axis equal; axis tight; axis off; title('mean non-rigid corrected','fontsize',14,'fontweight','bold')
    subplot(2,3,4); plot(1:T,cY,1:T,cM1,1:T,cM2); legend('raw data','rigid','non-rigid'); title('correlation coefficients','fontsize',14,'fontweight','bold')
    subplot(2,3,5); scatter(cY,cM1); hold on; plot([0.9*min(cY),1.05*max(cM1)],[0.9*min(cY),1.05*max(cM1)],'--r'); axis square;
        xlabel('raw data','fontsize',14,'fontweight','bold'); ylabel('rigid corrected','fontsize',14,'fontweight','bold');
    subplot(2,3,6); scatter(cM2,cM4); hold on; plot([0.9*min(cY),1.05*max(cM1)],[0.9*min(cY),1.05*max(cM1)],'--r'); axis square;
        xlabel('rigid corrected','fontsize',14,'fontweight','bold'); ylabel('non-rigid corrected','fontsize',14,'fontweight','bold');

%% plot shifts        

shifts_r = horzcat(shifts1(:).shifts)';
shifts_nr = cat(ndims(shifts2(1).shifts)+1,shifts2(:).shifts);
shifts_nr = reshape(shifts_nr,[],ndims(Y)-1,T);
shifts_x = squeeze(shifts_nr(:,1,:))';
shifts_y = squeeze(shifts_nr(:,2,:))';

patch_id = 1:size(shifts_x,2);
str = strtrim(cellstr(int2str(patch_id.')));
str = cellfun(@(x) ['patch # ',x],str,'un',0);



%%

figure;
for t = 1:1:T
    %subplot(131);imagesc(Y(:,:,t),[nnY,mmY]); xlabel('raw data','fontsize',14,'fontweight','bold'); axis equal; axis tight; axis off;  
    %subplot(132);imagesc(M1(:,:,t),[nnY,mmY]); xlabel('rigid corrected','fontsize',14,'fontweight','bold'); axis equal; axis tight;
    % set(gca,'XTick',[],'YTick',[]);
    %title(sprintf('Frame %i out of %i',t,T),'fontweight','bold','fontsize',14); colormap('bone')
    subplot(121);imagesc(Y(:,:,t),[nnY,mmY]); xlabel('non-rigid corrected','fontsize',14,'fontweight','bold'); axis equal; axis tight;
    title(sprintf('Frame %i out of %i',t,T),'fontweight','bold','fontsize',14); colormap('bone')
    subplot(122);imagesc(M1(:,:,t),[nnY,mmY]); xlabel('non-rigid corrected','fontsize',14,'fontweight','bold'); axis equal; axis tight;
    title(sprintf('Frame %i out of %i',t,T),'fontweight','bold','fontsize',14); colormap('bone')
    set(gca,'XTick',[],'YTick',[]);
    drawnow;
    pause(0.02);
end