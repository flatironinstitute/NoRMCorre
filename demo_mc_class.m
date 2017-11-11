% same as demo.m and demo_1p.m but using the MotionCorrection object. Scroll at the end
% to the 1p example

clear
gcp;

name = 'granule_love2.tif';         % two photon dataset
if ~exist(name,'file')  % download file if it doesn't exist in the directory
    url = 'https://www.dropbox.com/s/mjmtwn4pdgydkny/granule_love2.tif.zip?dl=1';
    filename = 'granule_love2.tif.zip';
    fprintf('downloading the file...');
    outfilename = websave(filename,url);
    fprintf('...done. Now unzipping...')
    unzip(filename);
    fprintf('done.');
end

%% rigid motion correction 

MC_rigid = MotionCorrection(name);
options_rigid = NoRMCorreSetParms('d1',MC_rigid.dims(1),'d2',MC_rigid.dims(2),'bin_width',50,'max_shift',15,'us_fac',50,'init_batch',200);
MC_rigid.motionCorrectSerial(options_rigid);  % can also try parallel

%% pw-rigid motion correction (in parallel)

MC_nonrigid = MotionCorrection(name);
options_nonrigid = NoRMCorreSetParms('d1',MC_nonrigid.dims(1),'d2',MC_nonrigid.dims(2),'grid_size',[32,32],'mot_uf',4,'bin_width',50,'max_shift',15,'max_dev',3,'us_fac',50,'init_batch',200);
MC_nonrigid.motionCorrectSerial(options_nonrigid);

%% compute metrics

MC_rigid.correlationMean(10);
MC_rigid.crispness(10);
MC_nonrigid.correlationMean(10);
MC_nonrigid.crispness(10);

%% do some plotting
nnY = quantile(MC_rigid.meanY(:),0.005);
mmY = quantile(MC_rigid.meanY(:),0.995);
T = MC_rigid.T;
cY = MC_rigid.corrY;
cM1 = MC_rigid.corrM;
cM2 = MC_nonrigid.corrM;
figure;
    ax1 = subplot(2,3,1); imagesc(MC_rigid.meanY,[nnY,mmY]);  axis equal; axis tight; axis off; title('mean raw data','fontsize',14,'fontweight','bold')
    ax2 = subplot(2,3,2); imagesc(MC_rigid.meanM,[nnY,mmY]);  axis equal; axis tight; axis off; title('mean rigid corrected','fontsize',14,'fontweight','bold')
    ax3 = subplot(2,3,3); imagesc(MC_nonrigid.meanM,[nnY,mmY]); axis equal; axis tight; axis off; title('mean non-rigid corrected','fontsize',14,'fontweight','bold')
    subplot(2,3,4); plot(1:T,cY,1:T,cM1,1:T,cM2); legend('raw data','rigid','non-rigid'); title('correlation coefficients','fontsize',14,'fontweight','bold')
    subplot(2,3,5); scatter(cY,cM1); hold on; plot([0.9*min(cY),1.05*max(cM1)],[0.9*min(cY),1.05*max(cM1)],'--r'); axis square;
        xlabel('raw data','fontsize',14,'fontweight','bold'); ylabel('rigid corrected','fontsize',14,'fontweight','bold');
    subplot(2,3,6); scatter(cM1,cM2); hold on; plot([0.9*min(cY),1.05*max(cM1)],[0.9*min(cY),1.05*max(cM1)],'--r'); axis square;
        xlabel('rigid corrected','fontsize',14,'fontweight','bold'); ylabel('non-rigid corrected','fontsize',14,'fontweight','bold');
    linkaxes([ax1,ax2,ax3],'xy')
    
%% plot shifts

patch_id = 1:size(MC_nonrigid.shifts_x,2);
str = strtrim(cellstr(int2str(patch_id.')));
str = cellfun(@(x) ['patch # ',x],str,'un',0);

figure;
    ax1 = subplot(311); plot(1:T,cY,1:T,cM1,1:T,cM2); legend('raw data','rigid','non-rigid'); title('correlation coefficients','fontsize',14,'fontweight','bold')
            set(gca,'Xtick',[])
    ax2 = subplot(312); plot(MC_nonrigid.shifts_x); hold on; plot(MC_rigid.shifts_x,'--k','linewidth',2); title('displacements along x','fontsize',14,'fontweight','bold')
            set(gca,'Xtick',[])
    ax3 = subplot(313); plot(MC_nonrigid.shifts_y); hold on; plot(MC_rigid.shifts_y,'--k','linewidth',2); title('displacements along y','fontsize',14,'fontweight','bold')
            xlabel('timestep','fontsize',14,'fontweight','bold')
    linkaxes([ax1,ax2,ax3],'x'); legend(str)


%%
%%%%%% 1p example following demo_1p.m  %%%%%%%%

clear;
gcp;

name = '/Users/epnevmatikakis/Documents/Ca_datasets/Miniscope/msCam13.avi';
MC_1p = MotionCorrection(name);
MC_1p.HPF(7,17);                                        % high pass filter 
options_1p = NoRMCorreSetParms('d1',MC_1p.dims(1),'d2',MC_1p.dims(2),'bin_width',50,'max_shift',20,'iter',1,'correct_bidir',false);
MC_1p.motionCorrectParallel(options_1p);                % correct the high pass filtered file
MC_1p.M = MC_1p.applyShifts(MC_1p.file_orig);           % apply shifts to the original file

%% plot mean images
figure;subplot(121); imagesc(mean(MC_1p.file_orig,3)); title('Mean image (raw)'); axis tight; axis equal;
       subplot(122); imagesc(mean(MC_1p.M,3)); title('Mean image (corrected)'); axis tight; axis equal;