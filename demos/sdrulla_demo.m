clear
name = '/Users/epnevmatikakis/Documents/Ca_datasets/Sdrulla/ExampleStack1.tif';
Y = bigread2(name); % read the file (optional, you can also pass the path in the function instead of Y)
Y = double(Y);      % convert to double precision 

%% set parameters (first try out rigid motion correction)

options.grid_size = [size(Y,1),size(Y,2)];  % size of patch in each direction
options.bin_width = 30;                     % number of bins after which you update template
options.mot_uf = 1;                         % upsampling factor for smaller patches
options.us_fac = 6;                         % upsampling factor for subpixel registration
options.method = {'median','mean'};         % averaging method for computing and updating templates
options.overlap_pre = 16;                   % amount of overlap for each patch
options.overlap_post = 16;                  % amount of overlap for each patch
options.iter = 1;                           % number of passes (set to 1 for one-pass online processing)
options.plot_flag = true;                  % flag for plotting results while correcting
options.memmap = false;                     % save output in a .mat file
options.iter = 1;
options.max_shift = 10;

%% perform motion correction
tic; [M1,shifts1,template1] = online_motion_correction_patches(Y,options); toc

%% now try non-rigid
options.grid_size = [128,128];                % size of patch in each direction
options.mot_uf = 8;                         % further upsample by a given factor
options.plot_flag = true;
options.iter = 1;
options.method = {'median','mean'}; 
tic; [M5,shifts5,template5] = online_motion_correction_patches(Y,options); toc

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
T = length(cY);
%%
figure;
    subplot(2,3,1); imagesc(mY,[nnY,mmY]);  axis equal; axis tight; axis off; title('mean raw data','fontsize',14,'fontweight','bold')
    subplot(2,3,2); imagesc(mM2,[nnY,mmY]);  axis equal; axis tight; axis off; title('mean rigid corrected','fontsize',14,'fontweight','bold')
    subplot(2,3,3); imagesc(mM3,[nnY,mmY]); axis equal; axis tight; axis off; title('mean non-rigid corrected','fontsize',14,'fontweight','bold')
    subplot(2,3,4); plot(1:T,cY,1:T,cM2,1:T,cM3); legend('raw data','rigid','non-rigid'); title('correlation coefficients','fontsize',14,'fontweight','bold')
    subplot(2,3,5); scatter(cY,cM1); hold on; plot([0.9*min(cY),1.05*max(cM1)],[0.9*min(cY),1.05*max(cM1)],'--r'); axis square;
        xlabel('raw data','fontsize',14,'fontweight','bold'); ylabel('rigid corrected','fontsize',14,'fontweight','bold');
    subplot(2,3,6); scatter(cM2,cM3); hold on; plot([0.9*min(cY),1.05*max(cM1)],[0.9*min(cY),1.05*max(cM1)],'--r'); axis square;
        xlabel('rigid corrected','fontsize',14,'fontweight','bold'); ylabel('non-rigid corrected','fontsize',14,'fontweight','bold');
%% plot

figure;
for t = 1:2:T
    subplot(121);imagesc(Y(:,:,t),[nnY,mmY]); xlabel('raw data','fontsize',14,'fontweight','bold'); axis equal; axis tight; axis off;  
    %subplot(132);imagesc(M1(:,:,t),[nnY,mmY]); xlabel('rigid corrected','fontsize',14,'fontweight','bold'); axis equal; axis tight;
    % set(gca,'XTick',[],'YTick',[]);
    %title(sprintf('Frame %i out of %i',t,T),'fontweight','bold','fontsize',14); colormap('bone')
    subplot(122);imagesc(M2(:,:,t),[nnY,mmY]); xlabel('non-rigid corrected','fontsize',14,'fontweight','bold'); axis equal; axis tight;
    title(sprintf('Frame %i out of %i',t,T),'fontweight','bold','fontsize',14); colormap('bone')
    set(gca,'XTick',[],'YTick',[]);
    drawnow;
    pause(0.02);
end

%%
%%
figure;
    subplot(2,3,1); imagesc(dY,[nnY,mmY]);  axis equal; axis tight; axis off; title('mean raw data','fontsize',14,'fontweight','bold')
    subplot(2,3,2); imagesc(dM1,[nnY,mmY]);  axis equal; axis tight; axis off; title('mean rigid corrected','fontsize',14,'fontweight','bold')
    subplot(2,3,3); imagesc(dM2,[nnY,mmY]); axis equal; axis tight; axis off; title('mean non-rigid corrected','fontsize',14,'fontweight','bold')
    subplot(2,3,4); plot(1:T,cdY,1:T,cdM1,1:T,cdM2); legend('raw data','rigid','non-rigid'); title('correlation coefficients','fontsize',14,'fontweight','bold')
    subplot(2,3,5); scatter(cdY,cdM1); hold on; plot([0.9*min(cY),1.05*max(cM1)],[0.9*min(cY),1.05*max(cM1)],'--r'); axis square;
        xlabel('raw data','fontsize',14,'fontweight','bold'); ylabel('rigid corrected','fontsize',14,'fontweight','bold');
    subplot(2,3,6); scatter(cdM1,cdM2); hold on; plot([0.9*min(cY),1.05*max(cM1)],[0.9*min(cY),1.05*max(cM1)],'--r'); axis square;
        xlabel('rigid corrected','fontsize',14,'fontweight','bold'); ylabel('non-rigid corrected','fontsize',14,'fontweight','bold');