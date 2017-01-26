%% taxidis demo
clear;
%name = '/Users/epnevmatikakis/Documents/Ca_datasets/Taxidis/mouse2/20160919_18_59_06__XYT.raw';
%name = '/Users/epnevmatikakis/Documents/Ca_datasets/Taxidis/single-channel recording/mouse2/20160919_18_59_06__XYT.raw';
%name = '/Users/epnevmatikakis/Documents/Ca_datasets/Taxidis/single-channel recording/mouse1/20160919_17_34_28__XYT.raw';
name = '/mnt/xfs1/home/eftychios/Documents/Ca_datasets/Taxidis/single-channel recording/mouse1/20160919_17_34_28__XYT.raw';
gcp;  % start parallel processing toolbox

%% perform rigid subpixel registration (runs in about 8 min on my machine)
options_rigid = NoRMCorreSetParms('d1',512,'d2',512,'bin_width',100,...
    'max_shift',32,'us_fac',50,'output_type','h5','h5_filename',[name(1:end-4),'.h5']);
tic; [M1,shifts1,template1] = normcorre_batch(name,options_rigid); toc

%% compare with other shifts
%load('/Users/epnevmatikakis/Documents/Ca_datasets/Taxidis/single-channel recording/mouse2/shifts.mat');
load('/Users/epnevmatikakis/Documents/Ca_datasets/Taxidis/single-channel recording/mouse1/shifts.mat');
T = length(shifts1);
shifts_tax = struct('shifts',cell(T,1),'shifts_up',cell(T,1),'diff',cell(T,1));
for t = 1:T
    shifts_tax(t).shifts = shifts(t,:)';
    shifts_tax(t).shifts_up = shifts(t,:)';
    shifts_tax(t).diff = 0;
end
options_tax = NoRMCorreSetParms('d1',512,'d2',512,'bin_width',100,...
    'max_shift',32,'us_fac',50,'output_type','h5','h5_filename',[name(1:end-4),'tax.h5']);
tic; M_tax = apply_shifts(name,shifts_tax,options_tax); toc;

%% perform non-rigid registration (much slower, runs in about 25 mins on my machine)
options_nonrigid = NoRMCorreSetParms('d1',512,'d2',512,'grid_size',[128,128],...
    'overlap_pre',64,'mot_uf',4,'bin_width',100,'max_shift',24,'max_dev',8,'us_fac',50,...
    'output_type','h5','h5_filename',[name(1:end-4),'_nr.h5']);
tic; [M2,shifts2,template2] = normcorre_batch(name,options_nonrigid); toc

%% load a chunk and compare the means
t_st = 1;           % first frame to read
num2read = 4000;    % number of frames to read
Y_raw = read_raw_file(name,t_st,num2read,[512,512],2);          % raw file
Y_tax = bigread2([name(1:end-4),'tax.h5'],t_st,num2read);       % shifts you gave me
Y_rig = bigread2([name(1:end-4),'.h5'],t_st,num2read);          % rigid correction
Y_nr = bigread2([name(1:end-4),'_nr.h5'],1,4000);               % non-rigid correction
%%
[c_raw,m_raw,v_raw] = motion_metrics(Y_raw,options_rigid.max_shift);
[c_tax,m_tax,v_tax] = motion_metrics(Y_tax,options_rigid.max_shift);
[c_rig,m_rig,v_rig] = motion_metrics(Y_rig,options_rigid.max_shift);
[c_nr, m_nr, v_nr ] = motion_metrics(Y_nr,options_rigid.max_shift);%figure;plot(1:4000,cY,1:4000,cM1,1:4000,cM2)

%% plot correlation coefficients
figure;plot(1:num2read,c_raw,1:num2read,c_tax,1:num2read,c_rig,1:num2read,c_nr); legend('raw','integer','subpixel','non-rigid');
        title('correlation coefficients');
%% plot mean and crispness values       
figure; colormap('bone')
    ax1 = subplot(221); imagesc(m_raw); axis square; axis off; title(sprintf('raw, crispness: %i',round(v_raw)))
    ax2 = subplot(222); imagesc(m_tax,[min(m_raw(:)),max(m_raw(:))]); axis square; axis off; title(sprintf('rigid integer, crispness: %i',round(v_tax)));
    ax3 = subplot(223); imagesc(m_rig,[min(m_raw(:)),max(m_raw(:))]); axis square; axis off; title(sprintf('rigid subpixel, crispness: %i',round(v_rig)));
    ax4 = subplot(224); imagesc(m_nr,[min(m_raw(:)),max(m_raw(:))]); axis square; axis off; title(sprintf('non-rigid, crispness: %i',round(v_nr)))
    linkaxes([ax1,ax2,ax3,ax4],'xy')

%% delete the h5 files     
%delete('/Users/epnevmatikakis/Documents/Ca_datasets/Taxidis/single-channel recording/mouse1/*.h5')