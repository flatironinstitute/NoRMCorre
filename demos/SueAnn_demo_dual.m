clear;
name = '/Users/epnevmatikakis/Documents/Ca_datasets/Sueann/Sue/20161121/k59_20161121_MMP_150um_47mW940nm_zoom4p2_00001_00029.tif';
name_gr = '/Users/epnevmatikakis/Documents/Ca_datasets/Sueann/Sue/k59_20161209_MMP_150um_47mW940nm_zoom2p2_00001_00002_g.tif';
name_rd = '/Users/epnevmatikakis/Documents/Ca_datasets/Sueann/Sue/k59_20161209_MMP_150um_47mW940nm_zoom2p2_00001_00002_r.tif';

addpath(genpath('../../NoRMCorre'));
Yg = read_file(name_gr);
Yr = read_file(name_rd);
%Y = read_file(name);
%Yr = Y(:,:,2:2:end);
%Yg = Y(:,:,1:2:end);
%%
% figure;
% for t = 1:3000
%     subplot(121); imagesc(Yr(:,:,t),[-300,-100]); axis square;
%     subplot(122); imagesc(Yg(:,:,t),[-400,2400]); axis square;
%     title(num2str(t))
%     drawnow;
%     pause(0.01)
% end

%% set parameters (first try out rigid motion correction)
[d1,d2,T] = size(Yr);
Yr = double(Yr);
Yg = double(Yg);
%%
options_r = NoRMCorreSetParms('d1',d1,'d2',d2,'bin_width',30,'max_shift',15,'iter',2);

tic; [M1_red,shifts1_red,template1_red] = normcorre_batch(Yr,options_r); toc % register filtered data
tic; [M1_grn,shifts1_grn,template1_grn] = normcorre_batch(Yg,options_r); toc % register filtered data

%%
[cYr,mYr,vYr] = motion_metrics(Yr,options_r.max_shift);
[cYg,mYg,vYg] = motion_metrics(Yg,options_r.max_shift);
[cM1r,mM1r,vM1r] = motion_metrics(M1_red,options_r.max_shift);
[cM1g,mM1g,vM1g] = motion_metrics(M1_grn,options_r.max_shift);
%%

figure;
    ax1 = subplot(121); imagesc(mM1r); axis square;
    ax2 = subplot(122); imagesc(mM1g); axis square;
    linkaxes([ax1,ax2],'xy')

%%    
shifts_red = horzcat(shifts1_red(:).shifts)';
shifts_grn = horzcat(shifts1_grn(:).shifts)';

%% now try non-rigid
options_nr = NoRMCorreSetParms('d1',d1,'d2',d2,'bin_width',30, ...
    'grid_size',[128,128],'mot_uf',4, ...
    'overlap_pre',32,'overlap_post',32,'max_shift',20);

tic; [M2_red,shifts2_red,template2_red] = normcorre_batch(Yr,options_nr); toc
tic; [M2_grn,shifts2_grn,template2_grn] = normcorre_batch(Yg,options_nr); toc % register filtered data

%%
tic; Mrg = apply_shifts(Yg,shifts2_red,options_nr); toc % apply the shifts to the removed percentile

[cM2rg,mM2rg,vM2rg] = motion_metrics(Mrg,options_nr.max_shift);

% These results are with only one pass of the data completely online. In
% practice we can improve if we do multiple passes (with options.iter > 1)
% or do fancy things like, first correct with rigid which is more robust to
% large perturbation, then refine by correcting with non-rigid

%%

[cM2r,mM2r,vM2r] = motion_metrics(M2_red,options_r.max_shift);
[cM2g,mM2g,vM2g] = motion_metrics(M2_grn,options_r.max_shift);

%%

figure;
    ax1 = subplot(121); imagesc(mM1g); axis square;
    ax2 = subplot(122); imagesc(mM2g); axis square;
    linkaxes([ax1,ax2],'xy')


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

shifts_r = horzcat(shifts1_grn(:).shifts)';
shifts_nr = cat(ndims(shifts2_grn(1).shifts)+1,shifts2_grn(:).shifts);
shifts_nr = reshape(shifts_nr,[],ndims(Yr)-1,T);
shifts_ng = cat(ndims(shifts2_red(1).shifts)+1,shifts2_red(:).shifts);
shifts_ng = reshape(shifts_ng,[],ndims(Yr)-1,T);

shifts_x = squeeze(shifts_nr(:,1,:))';
shifts_y = squeeze(shifts_nr(:,2,:))';

patch_id = 1:size(shifts_x,2);
str = strtrim(cellstr(int2str(patch_id.')));
str = cellfun(@(x) ['patch # ',x],str,'un',0);

figure;
    ax1 = subplot(311); plot(1:T,cYg,1:T,cM1g,1:T,cM2g); legend('raw data','rigid','non-rigid'); title('correlation coefficients','fontsize',14,'fontweight','bold')
            set(gca,'Xtick',[])
    ax2 = subplot(312); plot(shifts_x); hold on; plot(shifts_r(:,1),'--k','linewidth',2); title('displacements along x','fontsize',14,'fontweight','bold')
            set(gca,'Xtick',[])
    ax3 = subplot(313); plot(shifts_y); hold on; plot(shifts_r(:,2),'--k','linewidth',2); title('displacements along y','fontsize',14,'fontweight','bold')
            xlabel('timestep','fontsize',14,'fontweight','bold')
    linkaxes([ax1,ax2,ax3],'x')

%% display downsampled data

tsub = 10;

Y_ds = downsample_data(Yg,'time',tsub);
M_ds = downsample_data(M2_grn,'time',tsub);
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