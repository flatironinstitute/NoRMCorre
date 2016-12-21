clear;
name = '/Users/epnevmatikakis/Downloads/m11273A_00001_00001.tif';
addpath(genpath('../../NoRMCorre'));
tiffInfo = imfinfo(name);
T = length(tiffInfo);
Y1 = imread(name,'Index',1,'Info',tiffInfo);
Y = zeros([size(Y1),T]);
Y(:,:,1) = double(Y1);
tic;
for t = 2:T
    Y(:,:,t) = double(imread(name,'Index',t,'Info',tiffInfo));
end
toc

%% set parameters (first try out rigid motion correction)

options_r.grid_size = [size(Y,1),size(Y,2)];  % size of patch in each direction
options_r.bin_width = 30;                     % number of bins after which you update template
options_r.mot_uf = 1;                         % upsampling factor for smaller patches
options_r.us_fac = 10;                         % upsampling factor for subpixel registration
options_r.method = {'median','mean'};           % averaging method for computing and updating templates
options_r.overlap_pre = 16;                   % amount of overlap for each patch
options_r.overlap_post = 16;                  % amount of overlap for each patch
options_r.plot_flag = false;                  % flag for plotting results while correcting
options_r.memmap = false;                     % save output in a .mat file
options_r.iter = 1;
options_r.max_shift = 15;

%% perform motion correction
%profile on
tic; [M1,shifts1,template1] = normcorre_batch(Y,options_r); toc
%profile off
%profile viewer

%% now try non-rigid

options_nr = options_r;
options_nr.grid_size = [256,256]/2;              % size of patch in each direction
options_nr.mot_uf = 4;                         % further upsample by a given factor
options_nr.plot_flag = false;                  % flag for plotting results while correcting
options_nr.make_avi = false;
options_nr.overlap_pre = 64;
options_nr.overlap_post = 32;
options_nr.max_dev = 4;
options_nr.bin_width = 50;
%%
tic; [M2,shifts2,template2] = normcorre(Y2,options_nr,template2); toc

%%
Y2 = loadtiff('/Users/epnevmatikakis/Downloads/m11273A_00001_00002.tif');
tic; [M2b,shifts2b,template2] = normcorre_batch(Y2,options_nr,template2); toc

% These results are with only one pass of the data completely online. In
% practice we can improve if we do multiple passes (with options.iter > 1)
% or do fancy things like, first correct with rigid which is more robust to
% large perturbation, then refine by correcting with non-rigid

%% plot data

nnY = quantile(Y(1:10000),0.005);
mmY = quantile(Y(1:10000),0.995);
%%
[cY,mY,vY] = motion_metrics(Y,options_nr.max_shift);
[cM1,mM1,vM1] = motion_metrics(M2b,options_nr.max_shift);
[cM2,mM2,vM2] = motion_metrics(M2,options_nr.max_shift);

T = length(cY);
%% plot metrics
figure;
    ax1 = subplot(2,3,1); imagesc(mY,[nnY,mmY]);  axis equal; axis tight; axis off; title('mean raw data','fontsize',14,'fontweight','bold')
    ax2 = subplot(2,3,2); imagesc(mM1,[nnY,mmY]); axis equal; axis tight; axis off; title('mean rigid corrected','fontsize',14,'fontweight','bold')
    ax3 = subplot(2,3,3); imagesc(mM2,[nnY,mmY]); axis equal; axis tight; axis off; title('mean non-rigid corrected','fontsize',14,'fontweight','bold')
    subplot(2,3,4); plot(1:T,cY,1:T,cM1,1:T,cM2); legend('raw data','rigid','non-rigid'); title('correlation coefficients','fontsize',14,'fontweight','bold')
    subplot(2,3,5); scatter(cY,cM1); hold on; plot([0.9*min(cY),1.05*max(cM1)],[0.9*min(cY),1.05*max(cM1)],'--r'); axis square;
        xlabel('raw data','fontsize',14,'fontweight','bold'); ylabel('rigid corrected','fontsize',14,'fontweight','bold');
    subplot(2,3,6); scatter(cM1,cM2); hold on; plot([0.9*min(cY),1.05*max(cM1)],[0.9*min(cY),1.05*max(cM1)],'--r'); axis square;
        xlabel('rigid corrected','fontsize',14,'fontweight','bold'); ylabel('non-rigid corrected','fontsize',14,'fontweight','bold');
    linkaxes([ax1,ax2,ax3],'xy')
%% plot shifts        

shifts_r = horzcat(shifts1(:).shifts)';
shifts_nr = cat(ndims(shifts2(1).shifts)+1,shifts2(:).shifts);
shifts_nr = reshape(shifts_nr,[],ndims(Y)-1,T);
shifts_x = squeeze(shifts_nr(:,1,:))';
shifts_y = squeeze(shifts_nr(:,2,:))';

patch_id = 1:size(shifts_x,2);
str = strtrim(cellstr(int2str(patch_id.')));
str = cellfun(@(x) ['patch # ',x],str,'un',0);

figure;
    ax1 = subplot(311); plot(1:T,cY,1:T,cM1,1:T,cM2); legend('raw data','rigid','non-rigid'); title('correlation coefficients','fontsize',14,'fontweight','bold')
            set(gca,'Xtick',[])
    ax2 = subplot(312); plot(shifts_x); hold on; plot(shifts_r(:,1),'--k','linewidth',2); title('displacements along x','fontsize',14,'fontweight','bold')
            set(gca,'Xtick',[])
    ax3 = subplot(313); plot(shifts_y); hold on; plot(shifts_r(:,2),'--k','linewidth',2); title('displacements along y','fontsize',14,'fontweight','bold')
            xlabel('timestep','fontsize',14,'fontweight','bold')
    linkaxes([ax1,ax2,ax3],'x')

%% display downsampled data

tsub = 10;

Y_ds = downsample_data(Y,'time',tsub);
M_ds = downsample_data(M2b,'time',tsub);
nnY_ds = quantile(Y_ds(:),0.005);
mmY_ds = quantile(Y_ds(:),0.995);
%%
figure;
for t = 1:1:size(Y_ds,3)
    subplot(121);imagesc(Y_ds(:,:,t),[nnY_ds,mmY_ds]); xlabel('Raw data (downsampled)','fontsize',14,'fontweight','bold'); axis equal; axis tight;
    title(sprintf('Frame %i out of %i',t,T),'fontweight','bold','fontsize',14); colormap('bone');
    set(gca,'XTick',[],'YTick',[]);
    subplot(122);imagesc(M_ds(:,:,t),[nnY_ds,mmY_ds]); xlabel('non-rigid corrected','fontsize',14,'fontweight','bold'); axis equal; axis tight;
    title(sprintf('Frame %i out of %i',t,T),'fontweight','bold','fontsize',14); colormap('bone')
    set(gca,'XTick',[],'YTick',[]);
    drawnow;
    pause; %(0.02);
end