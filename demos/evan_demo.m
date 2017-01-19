clear
%foldername = '/Users/epnevmatikakis/Dropbox (Simons Foundation)/Evan/fly1_runA_Oct/';
foldername = '/mnt/xfs1/home/eftychios/Dropbox (Simons Foundation)/Evan/fly1_runA_Oct/';
name = [foldername,'R_fly1_runA_Oct_1.tiff'];

Z = loadtiff(name);
T = 300;        % number of volumes to read
Y = zeros([size(Z),T],'uint16');
Y(:,:,:,1) = Z;

%%

for t = 2:T
    name = [foldername,'R_fly1_runA_Oct_',num2str(t),'.tiff'];
    Z = loadtiff(name);
    Y(:,:,:,t) = Z;
    disp(t);
end

savefast('evan_fly1_runA_RC.mat.mat','Y');  %% save loaded file in a mat file
data = matfile('evan_fly1_runA_RC.mat.mat','Writable',true);
%%
sizY = size(data,'Y');
options = NoRMCorreSetParms('d1',sizY(1),'d2',sizY(2),'d3',sizY(3)','grid_size',[64,64,50],'bin_width',10,'mot_uf',2,'us_fac',10,...
            'method',{'median','mean'},'overlap_pre',16,'overlap_post',16,'max_dev',[4,4,4],...
            'output_type','h5','mem_batch_size',10,'h5_filename','evan_fly1_runA_RC_MC.h5','use_parallel',false,...
            'min_patch_size',[32,32,16]);
%options.name = [name,'_corrected.avi'];     % name for movie

%% perform motion correction
tic;
%profile on;
[M,shifts,template] = normcorre(data,options);
%profile off;
%profile viewer
toc
%%
% [cY,mY,vY] = motion_metrics(data,10);
% [cM,mM,vM] = motion_metrics(M,10);
% figure;plot(1:T,cY,1:T,cM)

%%
data_rc = matfile('evan_fly1_runA_RC_mc.mat','Writable',true);

%% load each volume and apply shifts
T = 300;
M = zeros(sizY,'single');
Y = M;
for t = 1:T
    name = ['/mnt/xfs1/home/eftychios/Dropbox (Simons Foundation)/Evan/fly1_runA_Oct_gc/G_fly1_runA_Oct_',num2str(t),'.tiff'];
    Z = loadtiff(name);
    Y(:,:,:,t) = Z;
    M(:,:,:,t) = apply_shifts(Z,shifts(t),options);
    disp(t);
end

%savefast('evan_fly1_runA_GC_mc.mat','Y','M','shifts');
%data = matfile('evan_fly1_runA_GC_mc.mat','Writable',true);

%% compute metrics (truncate some data for efficiency)
[cY,mY,vY] = motion_metrics(double(Y(50:130,100:430,40:90,:)),16);
[cM,mM,vM] = motion_metrics(double(M(50:130,100:430,40:90,:)),16);
figure;plot(1:T,cY,1:T,cM);

%% plot a plane

pl = 80;
Y_pl = squeeze(Y(pl,:,:,:));
M_pl = squeeze(M(pl,:,:,:));
nnY = 0; %quantile(Y(:),0.0025);
mmY = quantile(Y(:),0.9975);
%%

fig = figure;
    set(gcf, 'PaperUnits', 'points', 'Units', 'points');
    set(gcf, 'Position', round([100 100 2*sizY(2) 2*sizY(2)]));
while 1
for t = 1:1:T
    imagesc([Y_pl(:,:,t),M_pl(:,:,t)],[nnY,mmY]); axis equal; axis tight;
    ylabel('Y axis','fontweight','bold','fontsize',14);
    xlabel('Z axis','fontweight','bold','fontsize',14);
    title(sprintf('Frame %i out of %i',t,T),'fontweight','bold','fontsize',14); colormap('bone')    
    set(gca,'XTick',[],'YTick',[]);
    drawnow;
    pause(0.05);
end
end

%% save volumes
foldername = '/mnt/xfs1/home/eftychios/Dropbox (Simons Foundation)/Evan/fly1_runA_Oct_mc_RED/';
for t = 1:T
    filename = [foldername,'G_fly1_runA_Oct_',num2str(t),'_mc.tiff'];
    Z_temp = bigread2('evan_fly1_runA_RC_MC.h5',t,1);
    saveastiff(Z_temp, filename);
    disp(t)
end

%% save green volumes
mkdir('/mnt/xfs1/home/eftychios/Dropbox (Simons Foundation)/Evan','fly1_runA_Oct_mc_gc');
foldername_gc = '/mnt/xfs1/home/eftychios/Dropbox (Simons Foundation)/Evan/fly1_runA_Oct_mc_gc/';
for t = 1:T
    filename = [foldername_gc,'G_fly1_runA_Oct_',num2str(t),'_mc.tiff'];
    saveastiff(uint16(M(:,:,:,t)), filename);
    disp(t)
end