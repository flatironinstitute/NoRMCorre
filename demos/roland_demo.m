clear;
name1 = '/Users/epnevmatikakis/Documents/Ca_datasets/Zemla/243SP_dendrites_odorA.hdf5';
name2 = '/Users/epnevmatikakis/Documents/Ca_datasets/Zemla/243SP_dendrites_odorB.hdf5';

Y1 = read_file(name1);
Y2 = read_file(name2);

%% reshape data into 3d
[d1,d2,~] = size(Y1);
d3 = 6;
Y1 = reshape(Y1,d1,d2,d3,[]);
T1 = size(Y1,4);
Y2 = reshape(Y2,d1,d2,d3,[]);
T2 = size(Y2,4);
Y = cat(4,Y1,Y2);
T = size(Y,4);

%% motion correction
options = NoRMCorreSetParms('d1',d1,'d2',d2,'d3',d3,'grid_size',[d1,d2,d3],'bin_width',50,'mot_uf',[4,4,1],...
            'us_fac',50,'overlap_pre',16,'overlap_post',16);
% options.grid_size = [d1,d2,d3];              % size of patch in each direction
% options.bin_width = 50;                     % number of bins after which you update template
% options.mot_uf = 1;[4,4,1];                   % upsampling factor for smaller patches
% options.us_fac = 8;                         % upsampling factor for subpixel registration
% options.method = {'median','mean'};         % averaging method for computing and updating templates
% options.overlap_pre = 16;                   % amount of overlap for each patch
% options.overlap_post = 16;                  % amount of overlap for each patch
% options.iter = 1;                           % number of passes (set to 1 for one-pass online processing)
% options.plot_flag = false;                  % flag for plotting results while correcting
% options.memmap = false;                     % save output in a .mat file
% options.make_avi = false;                   % make a movie showing the results
% options.fr = 15;                            % frame rate for movie
% options.filename = 'roland1.mat';
% options.use_parallel = false;
%options.name = [name,'_corrected.avi'];     % name for movie

%% rigid motion correction
perm = randperm(T);
Y = Y(:,:,:,perm);
tt1 = tic;
[M1,shifts1,template1] = normcorre(Y,options);
toc(tt1)
M1(:,:,:,perm) = M1;
shifts1(perm)  = shifts1;

%% non rigid motion correction 

options_nr = NoRMCorreSetParms('d1',d1,'d2',d2,'d3',d3,'grid_size',[64,64,3],'bin_width',50,'mot_uf',[4,4,1],...
            'us_fac',50,'overlap_pre',32,'overlap_post',[16,16,2]);

tt1 = tic;
[M2,shifts2,template2] = normcorre_batch(Y,options_nr);
toc(tt1)
M2(:,:,:,perm) = M2;
Y(:,:,:,perm) = Y;
shifts2(perm)  = shifts2;

%% evaluate motion correction

[cY,mY,vY] = motion_metrics(Y,[10,10,10,10,0,0]);
[cM1,mM1,vM1] = motion_metrics(M1,[10,10,10,10,0,0]);
[cM2,mM2,vM2] = motion_metrics(M2,[10,10,10,10,0,0]);
figure;scatter(cY,cM1); hold on; plot([.5,.7],[.5,.7],'--r')
figure;scatter(cM1,cM2); hold on; plot([.5,.7],[.5,.7],'--r')

%%
M3 = permute(zeros(size(M2),'single'),[1,2,4,3]);
mM3 = cell(6,1);
for i = 1:6
    [M3(:,:,:,i),mM3{i}] = normcorre_batch(squeeze(Y(:,:,i,:)),options_nr);
end
%%
M3 = permute(M3,[1,2,4,3]);
[cM3,mM3,vM3] = motion_metrics(M3,[10,10,10,10,0,0]);
%% view a plane

pl = 1;randi(d3);
Ypl = squeeze(Y(:,:,pl,:));
Mpl = squeeze(M2(:,:,pl,:));
nnY = quantile(Ypl(:),0.01);
mmY = quantile(Ypl(:),0.99);
figure;
for t = 1:1:T   
    subplot(121);imagesc(Ypl(:,:,t),[nnY,mmY]); xlabel('Raw data','fontsize',14,'fontweight','bold'); axis equal; axis tight;
    title(sprintf('Plane %i',pl),'fontweight','bold','fontsize',14); colormap('bone')
    subplot(122);imagesc(Mpl(:,:,t),[nnY,mmY]); xlabel('Motion corrected','fontsize',14,'fontweight','bold'); axis equal; axis tight;
    title(sprintf('Frame %i out of %i',t,T),'fontweight','bold','fontsize',14); colormap('bone')
    set(gca,'XTick',[],'YTick',[]);
    drawnow;
    pause(0.02);
end

%% now run CNMF

addpath(genpath('/Users/epnevmatikakis/Documents/MATLAB/github/ca_source_extraction'));
Mf = M;
shifts_cat = cat(4,shifts(:).shifts);
max_disp = squeeze(max(max(max(shifts_cat,[],4),[],2),[],1));
min_disp = squeeze(min(min(min(shifts_cat,[],4),[],2),[],1));
M = Mf(round(max_disp(1)+1):round(d1+min_disp(1)),round(max_disp(2)+1):round(d2+min_disp(2)),:,:);
[d1,d2,d3,T] = size(M);
d = d1*d2*d3;
%% Set parameters

K = 200;                                          % number of components to be found
tau = [5,5,2];                                    % std of gaussian kernel (size of neuron) 
p = 0;                                            % order of autoregressive system (p = 0 no dynamics for slow imaging rate)
merge_thr = 0.95;                                 % merging threshold

options = CNMFSetParms(...                      
    'd1',d1,'d2',d2,'d3',d3,...                  % dimensions of datasets
    'search_method','dilate',...                 % search locations when updating spatial components
    'maxIter',15,...                             % number of NMF iterations during initialization
    'deconv_method','constrained_foopsi',...     % activity deconvolution method
    'temporal_iter',2,...                        % number of block-coordinate descent steps 
    'fudge_factor',0.98,...                      % bias correction for AR coefficients
    'merge_thr',merge_thr,...                    % merging threshold
    'gSig',tau,'nb',1 ...
    );
%% Data pre-processing

[P,M] = preprocess_data(M,p);
Cn = correlation_image_3D(M); % for large datasets change with reshape(P.sn,d1,d2,d3), %max(Y,[],3); %std(Y,[],3); % image statistic (only for display purposes)


%%

[Ain,Cin,bin,fin,center] = initialize_components(M,K,tau,options,P);  % initialize
ff = find(sum(Ain)<1e-3*mean(sum(Ain)));   % remove very small components
Ain(:,ff) = [];
Cin(ff,:) = [];
center(ff,:) = [];

%% display centers of found components
plotCenteroverY(Cn, center, [d1,d2,d3]);  % plot found centers against max-projections of background image


%% update spatial components
Mr = reshape(M,d,T);
%clear Y;
[A,b,Cin] = update_spatial_components(Mr,Cin,fin,Ain,P,options);

%% update temporal components
P.p = 0;
[C,f,P,S] = update_temporal_components(Mr,A,b,Cin,fin,P,options);
%[C_df,Df] = extract_DF_F(Mr,A,C,P,options);

%% plot components max projections
comp = randi(size(A,2));    % select a random component
atemp = reshape(full(A(:,comp)),d1,d2,d3);
figure;
    subplot(131); imagesc(max(atemp,[],3)); title('xy projection')
    subplot(132); plot(squeeze(max(atemp,[],1))); title('yz projection');
        xlabel(['component ',num2str(comp)]);
    subplot(133); plot(squeeze(max(atemp,[],2))); title('xz projection');

%% repeat (optional)
%[A2,b2,Cm] = update_spatial_components(Yr,Cm,f,Am,Pm,options);
%[C2,f2,P2,S2] = update_temporal_components(Yr,A2,b2,Cm,f,Pm,options);

plot_components_3D_GUI(M,A,C,b,f,Cn,options);

%% view a plane

M_rec = A*C + b*f;
%%
pl = 1;randi(d3);
%M_rec = reshape(M_rec,d1,d2,d3,T);
Ypl = squeeze(M_rec(:,:,pl,:));
Mpl = squeeze(M(:,:,pl,:));
nnY = quantile(Ypl(:),0.01);
mmY = quantile(Ypl(:),0.99);
figure;
for t = 1:1:T   
    subplot(121);imagesc(Ypl(:,:,t),[nnY,mmY]); xlabel('Raw data','fontsize',14,'fontweight','bold'); axis equal; axis tight;
    title(sprintf('Plane %i',pl),'fontweight','bold','fontsize',14); colormap('bone')
    subplot(122);imagesc(Mpl(:,:,t),[nnY,mmY]); xlabel('Motion corrected','fontsize',14,'fontweight','bold'); axis equal; axis tight;
    title(sprintf('Frame %i out of %i',t,T),'fontweight','bold','fontsize',14); colormap('bone')
    set(gca,'XTick',[],'YTick',[]);
    drawnow;
    pause(0.02);
end
