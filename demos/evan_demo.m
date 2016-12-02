clear
name = '/Users/epnevmatikakis/Dropbox (Simons Foundation)/Evan/fly1_runA_Oct/R_fly1_runA_Oct_1.tiff';

Z = bigread2(name);

Y = zeros([size(Z),300],'uint16');
Y(:,:,:,1) = Z;

%%

for t = 2:300
    name = ['/Users/epnevmatikakis/Dropbox (Simons Foundation)/Evan/fly1_runA_Oct/R_fly1_runA_Oct_',num2str(t),'.tiff'];
    Z = bigread2(name);
    Y(:,:,:,t) = Z;
    disp(t);
end

savefast('evan_fly1_runA_RC.mat.mat','Y');
data = matfile('evan_fly1_runA_RC.mat','Writable',true);
%%

options.grid_size = [64,64,50];             % size of patch in each direction
options.bin_width = 50;                     % number of bins after which you update template
options.mot_uf = 2;                         % upsampling factor for smaller patches
options.us_fac = 4;                         % upsampling factor for subpixel registration
options.method = {'median','mean'};         % averaging method for computing and updating templates
options.overlap_pre = 16;                   % amount of overlap for each patch
options.overlap_post = 16;
options.max_dev = [4,4,4];
options.iter = 1;                           % number of passes (set to 1 for one-pass online processing)
options.plot_flag = false;                  % flag for plotting results while correcting
options.memmap = true;                      % save output in a .mat file
options.make_avi = false;                   % make a movie showing the results
options.fr = 15;                            % frame rate for movie
options.filename = 'evan_fly1_runA_RC_mc.mat';
options.use_parallel = false;
options.min_patch_size = [32,32,16];
%options.name = [name,'_corrected.avi'];     % name for movie

%% perform motion correction
tic;
%profile on;
[M,shifts,template] = online_motion_correction_patches(data,options);
%profile off;
%profile viewer
toc
%%
[cY,mY,vY] = motion_metrics(data,10);
[cM,mM,vM] = motion_metrics(M,10);
figure;plot(1:300,cY,1:300,cM)

%%
data_rc = matfile('evan_fly1_runA_RC_mc.mat','Writable',true);

%% load and memory map green channel (pre-compute quantities)
overlap_post = options.overlap_post*ones(1,3);
[d1,d2,d3,T] = size(data_rc,'M');
[xx_s,xx_f,yy_s,yy_f,zz_s,zz_f,xx_us,xx_uf,yy_us,yy_uf,zz_us,zz_uf] = construct_grid(options.grid_size,options.mot_uf,d1,d2,d3,options.min_patch_size);
nd = 2 + (d3>1);

temp_cell = mat2cell_ov(data_rc.template,xx_us,xx_uf,yy_us,yy_uf,zz_us,zz_uf,overlap_post,[d1,d2,d3]);
Nr = cell([length(xx_us),length(yy_us),length(zz_us)]);
Nc = cell(size(Nr));
Np = cell(size(Nr));
Bs = cell(size(Nr));
for i = 1:length(xx_us)
    for j = 1:length(yy_us)
        for k = 1:length(zz_us)
            [nr,nc,np] = size(temp_cell{i,j,k});
            nr = ifftshift(-fix(nr/2):ceil(nr/2)-1);
            nc = ifftshift(-fix(nc/2):ceil(nc/2)-1);
            np = ifftshift(-fix(np/2):ceil(np/2)-1);
            [Nc{i,j,k},Nr{i,j,k},Np{i,j,k}] = meshgrid(nc,nr,np);
            extended_grid = [max(xx_us(i)-overlap_post(1),1),min(xx_uf(i)+overlap_post(1),d1),max(yy_us(j)-overlap_post(2),1),min(yy_uf(j)+overlap_post(2),d2),max(zz_us(k)-overlap_post(3),1),min(zz_uf(k)+overlap_post(3),d3)];            
            Bs{i,j,k} = permute(construct_weights([xx_us(i),xx_uf(i),yy_us(j),yy_uf(j),zz_us(k),zz_uf(k)],extended_grid),[2,1,3]); 
        end
    end
end
if nd == 2; Np = cellfun(@(x) 0,Nr,'un',0); end
%if nd == 3; Bs = cellfun(@(x) ones(size(x)), Nr, 'un',0); end

%% 
name = '/Users/epnevmatikakis/Dropbox (Simons Foundation)/Evan/fly1_runA_Oct_gc/G_fly1_runA_Oct_1.tiff';

Z = bigread2(name);

Y = zeros([size(Z),300],'uint16');
M = Y;
Y(:,:,:,1) = Z;
M_temp = uint16(apply_shifts(double(Z),shifts(1).shifts_up,overlap_post,xx_us,xx_uf,yy_us,yy_uf,zz_us,zz_uf,options.us_fac,Nr,Nc,Np,Bs));
data.M(:,:,:,1) = M_temp;
%%

for t = 2:300
    name = ['/Users/epnevmatikakis/Dropbox (Simons Foundation)/Evan/fly1_runA_Oct_gc/G_fly1_runA_Oct_',num2str(t),'.tiff'];
    Z = bigread2(name);
    Y(:,:,:,t) = Z;
    data.M(:,:,:,t) = uint16(apply_shifts(double(Z),shifts(t).shifts_up,overlap_post,xx_us,xx_uf,yy_us,yy_uf,zz_us,zz_uf,options.us_fac,Nr,Nc,Np,Bs));
    disp(t);
end

%savefast('evan_fly1_runA_GC_mc.mat','Y','M','shifts');
%data = matfile('evan_fly1_runA_GC_mc.mat','Writable',true);

%%
[cY,mY,vY] = motion_metrics(data,10);
[cM,mM,vM] = motion_metrics(data,10);

%%

pl = 40;
Y = zeros(186,532,300);
M = Y;
for t = 1:300
    Y(:,:,t) = double(data.Y(:,:,pl,t));
    M(:,:,t) = double(data.M(:,:,pl,t));
end

%%

nnY = quantile(Y(:),0.025);
mmY = quantile(Y(:),0.975);
%%
figure;
for t = 1:1:T
    subplot(211);imagesc(Y(:,:,t),[nnY,mmY]); xlabel('raw data','fontsize',14,'fontweight','bold'); axis equal; axis tight;
    title(sprintf('Frame %i out of %i',t,T),'fontweight','bold','fontsize',14); colormap('bone')
    subplot(212);imagesc(M(:,:,t),[nnY,mmY]); xlabel('non-rigid corrected','fontsize',14,'fontweight','bold'); axis equal; axis tight;
    title(sprintf('Frame %i out of %i',t,T),'fontweight','bold','fontsize',14); colormap('bone')
    set(gca,'XTick',[],'YTick',[]);
    drawnow;
    pause(0.04);
end
    