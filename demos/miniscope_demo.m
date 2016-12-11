clear;
%% read data and convert to double
name = '/Users/epnevmatikakis/Downloads/H16_M17_S54/msCam5.avi';
addpath(genpath('../../NoRMCorre'));
Y = read_file(name);
Y = double(Y);
%% remove a percentile value
p = 10;
Yf = Y;
Ypc = prctile(Y,p,3);
Y = bsxfun(@minus,Yf,Ypc);
%Y = prctfilt(Yf,p,200,100);
%Ypc = bsxfun(@minus,Yf,Y);

%% first try out rigid motion correction
[d1,d2,T] = size(Y);
options_r = NoRMCorreSetParms('d1',d1,'d2',d2,'bin_width',30,'max_shift',20,'iter',2);

%% register data and apply shifts to removed percentile
tic; [M1,shifts1,template1,options_r] = normcorre(Y,options_r); toc % register filtered data
tic; I = apply_shifts(Ypc,shifts1,options_r); toc % apply shifts to removed percentile
Mr = M1 + I; % rigid motion

%% compute metrics
[cYa,mYa,vYa] = motion_metrics(Y,options_r.max_shift);
[cM1a,mM1a,vM1a] = motion_metrics(M1,options_r.max_shift);
[cYb,mYb,vYb] = motion_metrics(bsxfun(@plus,Y,Ypc),options_r.max_shift);
[cM1b,mM1b,vM1b] = motion_metrics(Mr,options_r.max_shift);

%% now apply non-rigid motion correction
options_nr = NoRMCorreSetParms('d1',d1,'d2',d2,'bin_width',30, ...
    'grid_size',[128,128],'mot_uf',4, ...
    'overlap_pre',32,'overlap_post',32,'max_shift',20);

tic; [M2,shifts2,template2,options_nr] = normcorre(Y,options_nr); toc % register filtered data
tic; I2 = apply_shifts(Ypc,shifts2,options_nr); toc % apply the shifts to the removed percentile
Mpr = M2 + I2;

%% compute metrics
[cM2a,mM2a,vM2a] = motion_metrics(M2,options_nr.max_shift);
[cM2b,mM2b,vM2b] = motion_metrics(Mpr,options_nr.max_shift);

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
    ax1 = subplot(311); plot(1:T,cYa,1:T,cM1a,1:T,cM2a); legend('raw data','rigid','non-rigid'); title('correlation coefficients','fontsize',14,'fontweight','bold')
            set(gca,'Xtick',[],'XLim',[0,T-3])
    ax2 = subplot(312); plot(shifts_x); hold on; plot(shifts_r(:,1),'--k','linewidth',2); title('displacements along x','fontsize',14,'fontweight','bold')
            set(gca,'Xtick',[])
    ax3 = subplot(313); plot(shifts_y); hold on; plot(shifts_r(:,2),'--k','linewidth',2); title('displacements along y','fontsize',14,'fontweight','bold')
            xlabel('timestep','fontsize',14,'fontweight','bold')
    linkaxes([ax1,ax2,ax3],'x')
    
%% display downsampled data

tsub = 5;

Y_ds = downsample_data(Y,'time',tsub);
Yf_ds = downsample_data(Yf,'time',tsub);
M1_ds = downsample_data(M1,'time',tsub);
M1f_ds = downsample_data(Mr,'time',tsub);
M2_ds = downsample_data(M2,'time',tsub);
M2f_ds = downsample_data(Mpr,'time',tsub);
nnY_ds = quantile(Y_ds(:),0.0005);
mmY_ds = quantile(Y_ds(:),0.9995);
nnYf_ds = quantile(Yf_ds(:),0.0005);
mmYf_ds = quantile(Yf_ds(:),0.99995);
%% 
make_avi = false;
if make_avi
    vidObj = VideoWriter('full.avi');
    set(vidObj,'FrameRate',30);
    open(vidObj);
end
fig = figure;
    screensize = get(0,'Screensize' );
    fac = min(min((screensize(3:4)-100)./[3*d2,d1]),10);
    set(gcf, 'PaperUnits', 'points', 'Units', 'points');
    set(gcf, 'Position', round([100 100 fac*3*d2 fac*d1]));

for t = 1:1:size(Y_ds,3)-1
    subplot(131);imagesc(Yf_ds(:,:,t),[nnYf_ds,mmYf_ds]); xlabel('Raw data (downsampled)','fontsize',14,'fontweight','bold'); axis equal; axis tight;
    %title(sprintf('Frame %i out of %i',t,size(Y_ds,3)),'fontweight','bold','fontsize',14); 
    colormap('bone');
    set(gca,'XTick',[],'YTick',[]);
    subplot(132);imagesc(M1f_ds(:,:,t),[nnYf_ds,mmYf_ds]); xlabel('rigid corrected','fontsize',14,'fontweight','bold'); axis equal; axis tight;
    title(sprintf('Frame %i out of %i',t,size(Y_ds,3)),'fontweight','bold','fontsize',14); 
    colormap('bone')
    set(gca,'XTick',[],'YTick',[]);
    subplot(133);imagesc(M2f_ds(:,:,t),[nnYf_ds,mmYf_ds]); xlabel('non-rigid corrected','fontsize',14,'fontweight','bold'); axis equal; axis tight;
    %title(sprintf('Frame %i out of %i',t,size(Y_ds,3)),'fontweight','bold','fontsize',14); 
    colormap('bone')
    set(gca,'XTick',[],'YTick',[]);
%     subplot(234);imagesc(Yf_ds(:,:,t),[nnYf_ds,mmYf_ds]); axis equal; axis tight;
%     colormap('bone');
%     set(gca,'XTick',[],'YTick',[]);
%     subplot(235);imagesc(M1f_ds(:,:,t),[nnYf_ds,mmYf_ds]); axis equal; axis tight;
%     xlabel(sprintf('Frame %i out of %i',t,size(Y_ds,3)),'fontweight','bold','fontsize',14); colormap('bone')
%     set(gca,'XTick',[],'YTick',[]);
%     subplot(236);imagesc(M2f_ds(:,:,t),[nnYf_ds,mmYf_ds]); axis equal; axis tight;
%     colormap('bone')
%     set(gca,'XTick',[],'YTick',[]);
    drawnow;
    if make_avi  
        currFrame = getframe(fig);
        writeVideo(vidObj,currFrame);    
    end
end
if make_avi
    close(vidObj);
end
