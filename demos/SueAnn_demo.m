clear;
gcp;
name = '/Users/epnevmatikakis/Documents/Ca_datasets/Sueann/k56_20160608_RSM_125um_41mW_zoom2p2_00001_00034.tif';
addpath(genpath('../../NoRMCorre'));
Y = read_file(name);

%% set parameters (first try out rigid motion correction)
Y = Y(:,:,1:2000);
options_rigid = NoRMCorreSetParms('d1',size(Y,1),'d2',size(Y,2),'bin_width',50,'max_shift',15,'us_fac',50);

%% perform motion correction
tic; [M1,shifts1,template1] = normcorre_batch(Y,options_rigid); toc
[cY,mY,vY] = motion_metrics(Y,options_rigid.max_shift);
[cM1,mM1,vM1] = motion_metrics(M1,options_rigid.max_shift);

%% now try non-rigid motion correction
options_nonrigid = NoRMCorreSetParms('d1',size(Y,1),'d2',size(Y,2),'grid_size',[64,64],'overlap_pre',32,'overlap_post',32,'mot_uf',4,'bin_width',100,'max_shift',15,'max_dev',8,'us_fac',50);

tic; [M2,shifts2,template2] = normcorre_batch(Y,options_nonrigid,template1); toc
[cM2,mM2,vM2] = motion_metrics(M2,12);

%% plot data

nnY = quantile(Y(1:1000000),0.002);
mmY = quantile(Y(1:1000000),0.998);
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

%% make the movies of the paper
% %% make video of downsampled data
% tsub = 5;
% Y_sub = downsample_data(Y,'time',5);
% Mr_sub = downsample_data(M1,'time',5);
% Mn_sub = downsample_data(M2,'time',5);
% 
% nnY = quantile(Y_sub(1:1e7),0.005);
% mmY = quantile(Y_sub(1:1e7),1-0.005);
% 
% Y_all = cat(1,Y_sub,int16(Mr_sub),int16(Mn_sub));
% 
% fig = figure; colormap('bone')
%     set(gcf, 'Position', round([100 100 .8*512, 3*512]));
%     vidObj = VideoWriter('SueAnn_downsampled25xT.avi');
%     set(vidObj,'FrameRate',15);
%     open(vidObj);
% for t = 1:size(Y_sub,3)
%     imagesc(Y_all(:,:,t),[nnY,mmY]); axis equal; axis tight;
%         set(gca,'Xtick',[],'Ytick',[]);
%         set(gca,'Position',[0.05,0.05,0.92,0.92])
%         xlabel(sprintf('time %2.2f s',t/6),'fontsize',18,'fontweight','bold')
%         ylabel('Registered (pw-rigid)             Registered (rigid)                            Original data','fontsize',18,'fontweight','bold')
%     drawnow;
%     currFrame = getframe(fig);
%     writeVideo(vidObj,currFrame);    
% end
% close(vidObj);
% 
% %% make movie of online correction
% 
% options_nonrigid = NoRMCorreSetParms('d1',size(Y,1),'d2',size(Y,2),'grid_size',[64,64],'overlap_pre',32,'overlap_post',32,'mot_uf',4,'bin_width',50,'max_shift',15,'max_dev',3,'us_fac',50,'plot_flag',true,...
%     'make_avi',true,'name','SueAnnFullMovie.avi','fr',30);
% 
% tic; [M2b,shifts2b,template2b] = normcorre(Y,options_nonrigid); toc