%% taxidis demo
clear;
name = '/Users/epnevmatikakis/Documents/Ca_datasets/Taxidis/mouse2/20160919_18_59_06__XYT.raw';
gcp;
options_rigid = NoRMCorreSetParms('d1',512,'d2',512,'bin_width',100,...
    'max_shift',32,'us_fac',50,'output_type','h5','h5_filename',[name(1:end-4),'.h5']);
tic; [M1,shifts1,template1] = normcorre_batch(name,options_rigid); toc
%% apply other shifts
load('/Users/epnevmatikakis/Documents/Ca_datasets/Taxidis/mouse2/shifts_tax.mat');
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
%%

%% non-rigid
options_nonrigid = NoRMCorreSetParms('d1',512,'d2',512,'grid_size',[128,128],...
    'overlap_pre',64,'mot_uf',4,'bin_width',100,'max_shift',24,'max_dev',8,'us_fac',50,...
    'output_type','h5','h5_filename',[name(1:end-4),'_nr.h5']);
tic; [M2,shifts2,template2] = normcorre_batch(name,options_nonrigid); toc
%%
Y_raw = bigread2([name(1:end-4),'tax.h5'],1,4000);
Y_rig = bigread2([name(1:end-4),'.h5'],1,4000);
Y_nr = bigread2([name(1:end-4),'_nr.h5'],1,4000);
Y_nr2 = read_raw_file([name(1:end-4),'_nr.raw'],1,4000,[512,512],2);
[cY,mY,vY] = motion_metrics(Y_raw,options_rigid.max_shift);
[cM1,mM1,vM1] = motion_metrics(Y_rig,options_rigid.max_shift);
[cM2,mM2,vM2] = motion_metrics(Y_nr,options_rigid.max_shift);
figure;plot(1:4000,cY,1:4000,cM1,1:4000,cM2)
%%
figure;
    ax1 = subplot(131); imagesc(mY); axis square;
    ax2 = subplot(132); imagesc(mM1,[min(mY(:)),max(mY(:))]); axis square;
    ax3 = subplot(133); imagesc(mM2,[min(mY(:)),max(mY(:))]); axis square;
    linkaxes([ax1,ax2,ax3],'xy')
    
%%
name_in = [name(1:end-4),'_nr.h5'];
name_out = [name_in(1:end-2),'bin'];
h5_2_bin(name_in,name_out,'uint16',1000);