function shifts = register_frame_gpu(fftYMat,fftYPatches,fftTempMat,fftTempPatches,patches,options)

nd = ndims(fftYMat);
n_patches = [length(unique(patches(:,1))),length(unique(patches(:,3))),length(unique(patches(:,5)))];
shifts = zeros([n_patches,nd],'gpuArray');
N_patches = size(patches,1);

if N_patches > 1      
    if nd == 2; out_rig = dftregistration_min_max_gpu(fftTempMat,fftYMat,options.us_fac,-options.max_shift,options.max_shift,options.phase_flag); lb = out_rig(3:4); ub = out_rig(3:4); end
    if nd == 3; out_rig = dftregistration_min_max_3d(fftTempMat,fftYMat,options.us_fac,-options.max_shift,options.max_shift,options.phase_flag); lb = out_rig(3:5); ub = out_rig(3:5); end
    max_dev = options.max_dev;
else
    lb = -options.max_shift(1,nd);
    ub = options.max_shift(1,nd);
    max_dev = 0*options.max_dev;
end

for i = 1:N_patches
    [ix,iy,iz] = ind2sub(n_patches,i);
    if nd == 2
        output = dftregistration_min_max_gpu(fftTempPatches(:,:,:,i),fftYPatches(:,:,:,i),options.us_fac,lb-max_dev(1:2),ub+max_dev(1:2),options.phase_flag);  
    elseif nd == 3
        output = dftregistration_min_max_3d(fftTempPatches(:,:,:,i),fftYPatches(:,:,:,i),options.us_fac,lb-max_dev,ub+max_dev,options.phase_flag); 
        shifts(ix,iy,iz,3) = output(5);
    end

    shifts(ix,iy,iz,1) = output(3);
    shifts(ix,iy,iz,2) = output(4); 
end

end