clear;
name = '/Users/epnevmatikakis/Documents/Ca_datasets/Sueann/k56_20160608_RSM_125um_41mW_zoom2p2_00001_00034.tif';
addpath(genpath('../../NoRMCorre'));
Y = read_file(name);

%% set parameters (first try out rigid motion correction)

options_rigid = NoRMCorreSetParms('d1',size(Y,1),'d2',size(Y,2),'bin_width',50,'max_shift',15,'us_fac',50);

%% perform motion correction
tic; [M1a,shifts1a,template1a] = normcorre(Y,options_rigid); toc
tic; [M1b,shifts1b,template1b] = normcorre_batch(Y,options_rigid); toc
[cY,mY,vY] = motion_metrics(Y,options_rigid.max_shift);
[cM1a,mM1a,vM1a] = motion_metrics(M1a,options_rigid.max_shift);
[cM1b,mM1b,vM1b] = motion_metrics(M1b,options_rigid.max_shift);

%% now try non-rigid

%% now try non-rigid motion correction
options_nonrigid = NoRMCorreSetParms('d1',size(Y,1),'d2',size(Y,2),'grid_size',[128,128],'overlap_pre',32,'overlap_post',32,'mot_uf',4,'bin_width',50,'max_shift',15,'max_dev',3,'us_fac',50);

tic; [M2a,shifts2a,template2a] = normcorre(Y,options_nonrigid); toc
tic; [M2b,shifts2b,template2b] = normcorre_batch(Y,options_nonrigid); toc
[cM2a,mM2a,vM2a] = motion_metrics(M2a,options_rigid.max_shift);
[cM2b,mM2b,vM2b] = motion_metrics(M2b,options_rigid.max_shift);


% These results are with only one pass of the data completely online. In
% practice we can improve if we do multiple passes (with options.iter > 1)
% or do fancy things like, first correct with rigid which is more robust to
% large perturbation, then refine by correcting with non-rigid

%% plot data

nnY = quantile(Y(1:10000),0.005);
mmY = quantile(Y(1:10000),0.995);
%%
[cY,mY,vY] = motion_metrics(Y,options_nr.max_shift);
[cM1,mM1,vM1] = motion_metrics(M1,options_nr.max_shift);
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

tsub = 5;

Y_ds = downsample_data(Y,'time',tsub);
M_ds = downsample_data(M2,'time',tsub);
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
    pause(0.001);
end