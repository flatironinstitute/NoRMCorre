function [shifts,Mf] = register_frame(Yt,fftTempMat,fftTempPatches,patches,options)

minY = min(Yt(:));
maxY = max(Yt(:));
dim = size(Yt);
Yc = split_frame(Yt,patches);
fftY = fft(fft(fft(Yc,[],1),[],2),[],3);
nd = ndims(Yt);
n_patches = [length(unique(patches(:,1))),length(unique(patches(:,3))),length(unique(patches(:,5)))];
shifts = zeros([n_patches,nd]);
N_patches = size(patches,1);

if N_patches > 1      
    if nd == 2; out_rig = dftregistration_min_max(fftTempMat,fftn(Yt),options.us_fac,-options.max_shift,options.max_shift,options.phase_flag); lb = out_rig(3:4); ub = out_rig(3:4); end
    if nd == 3; out_rig = dftregistration_min_max_3d(fftTempMat,fftn(Yt),options.us_fac,-options.max_shift,options.max_shift,options.phase_flag); lb = out_rig(3:5); ub = out_rig(3:5); end
    max_dev = options.max_dev;
else
    lb = -options.max_shift(1,nd);
    ub = options.max_shift(1,nd);
    max_dev = 0*options.max_dev;
end

for i = 1:N_patches
    [ix,iy,iz] = ind2sub(n_patches,i);
    if nd == 2
        output = dftregistration_min_max(fftTempPatches(:,:,:,i),fftY(:,:,:,i),options.us_fac,lb-max_dev(1:2),ub+max_dev(1:2),options.phase_flag);  
    elseif nd == 3
        output = dftregistration_min_max_3d(fftTempPatches(:,:,:,i),fftY(:,:,:,i),options.us_fac,lb-max_dev,ub+max_dev,options.phase_flag); 
        shifts(ix,iy,iz,3) = output(5);
    end

    shifts(ix,iy,iz,1) = output(3);
    shifts(ix,iy,iz,2) = output(4); 
end

if nd == 3                
    shifts_up = zeros([dim,3]);
    do = size(shifts_up)./size(shifts);
    ds = size(shifts);
    if (0)
        [Xq,Yq,Zq] = meshgrid(linspace((2+0*1/do(2))/2,ds(2)+0*(0-1/do(2))/2,dim(2)),linspace((2+0*1/do(1))/2,ds(1)+0*(0-1/do(1))/2,dim(1)),linspace((2+0*1/do(3))/2,ds(3)+0*(0-1/do(3))/2,dim(3)));
    else
        [Xq,Yq,Zq] = meshgrid(linspace((1+1/do(2))/2,ds(2)+(1-1/do(2))/2,dim(2)),linspace((1+1/do(1))/2,ds(1)+(1-1/do(1))/2,dim(1)),linspace((1+1/do(3))/2,ds(3)+(1-1/do(3))/2,dim(3)));
        Xq(Xq<1)=1; Xq(Xq>dim(2))=dim(2);
        Yq(Yq<1)=1; Yq(Yq>dim(1))=dim(1);
        Zq(Zq<1)=1; Zq(Zq>dim(3))=dim(3);
    end
    if numel(shifts) > 3
%       tform = affine3d(diag([options.mot_uf(:);1]));
%       for dm = 1:3; shifts_up(:,:,:,dm) = imwarp(shifts(:,:,:,dm),tform,'OutputView',imref3d([options.d1,options.d2,options.d3])); end
        for dm = 1:3; shifts_up(:,:,:,dm) = interp3(shifts(:,:,:,dm),Xq,Yq,Zq,'makima'); end
    else
        for dm = 1:3; shifts_up(:,:,:,dm) = shifts(dm); end
    end
    shifts_up(2:2:end,:,:,2) = shifts_up(2:2:end,:,:,2) + options.col_shift;
    Mf = imwarp(Yt,-cat(4,shifts_up(:,:,:,2),shifts_up(:,:,:,1),shifts_up(:,:,:,3)),'bicubic','FillValues',options.add_value); 
else
    shifts_up = imresize(shifts,[options.d1,options.d2]);
    shifts_up(2:2:end,:,2) = shifts_up(2:2:end,:,2) + options.col_shift;
    Mf = imwarp(Yt,-cat(3,shifts_up(:,:,2),shifts_up(:,:,1)),'bicubic','FillValues',options.add_value);  
end   
Mf(Mf<minY)=minY;
Mf(Mf>maxY)=maxY;

end