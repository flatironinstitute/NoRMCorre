function [M,shifts,shifts_up] = register_frame(Yt,us_fac,template,grid_size,mot_uf,min_patch_size,max_shift,overlap_pre,overlap_post)

[d1,d2,d3] = size(Yt);
sizY = [d1,d2,d3];
if d3 == 1; nd = 2; else nd = 3; end
[xx_s,xx_f,yy_s,yy_f,zz_s,zz_f,xx_us,xx_uf,yy_us,yy_uf,zz_us,zz_uf] = construct_grid(grid_size,mot_uf,d1,d2,d3,min_patch_size);
Yc = mat2cell_ov(Yt,xx_s,xx_f,yy_s,yy_f,zz_s,zz_f,overlap_pre,sizY);
fftY = cellfun(@fftn, Yc, 'un',0);
temp_cell = mat2cell_ov(template,xx_s,xx_f,yy_s,yy_f,zz_s,zz_f,overlap_pre,sizY);
fftTemp = cellfun(@fftn,temp_cell,'un',0);

M_fin = cell(length(xx_us),length(yy_us),length(zz_us)); %zeros(size(Y_temp));
shifts_temp = zeros(length(xx_s),length(yy_s),length(zz_s),nd); 
for i = 1:length(xx_s)
    for j = 1:length(yy_s)           
        for k = 1:length(zz_s)
            if nd == 2
                [output,Greg] = dftregistration_max(fftTemp{i,j,k},fftY{i,j,k},us_fac,max_shift);                            
            elseif nd == 3
                [output,Greg] = dftregistration_max_3d(fftTemp{i,j,k},fftY{i,j,k},us_fac,max_shift);
                shifts_temp(i,j,k,3) = output(5);
            end
            shifts_temp(i,j,k,1) = output(3);
            shifts_temp(i,j,k,2) = output(4); 
            if mot_uf == 1
                M_temp = abs(ifftn(Greg));
                M_fin{i,j,k} = remove_boundaries(M_temp,output(3:end),'zero',template{i,j,k});
            end                                               
        end
    end
end            
   
shifts = squeeze(shifts_temp);

if any(mot_uf > 1)
    shifts_up = imresize(shifts_temp,[length(xx_uf),length(yy_uf)]);
    if mot_uf(3) > 1
        shifts_up = reshape(imresize(reshape(shifts_up,[length(xx_uf)*length(yy_uf),length(zz_f),nd]),[length(xx_uf)*length(yy_uf),length(zz_uf)]),[length(xx_uf),length(yy_uf),length(zz_uf),nd]);
    end
    for i = 1:length(xx_uf)
        for j = 1:length(yy_uf)
            for k = 1:length(zz_uf)
                extended_grid = [max(xx_us(i)-overlap_post,1),min(xx_uf(i)+overlap_post,d1),max(yy_us(j)-overlap_post,1),min(yy_uf(j)+overlap_post,d2),max(zz_us(k)-overlap_post,1),min(zz_uf(k)+overlap_post,d3)];
                I_temp = Yt(extended_grid(1):extended_grid(2),extended_grid(3):extended_grid(4),extended_grid(5):extended_grid(6));
                M_fin{i,j,k} = shift_reconstruct(I_temp,shifts_up(i,j,k,:),us_fac); %,Nr{i,j,k},Nc{i,j,k},Np{i,j,k});
            end
        end
    end
else
    shifts_up = shifts_temp;
end

%         ss = squeeze(shifts_up(:,:,t,:));
%         SS = mat2cell(ss,ones(size(ss,1),1),ones(size(ss,2),1),2);
%         Ic = mat2cell_ov(Yt,grid_size_fine,overlap_post);
%         M_fin = cellfun(@(x,s) shift_reconstruct(x,s,us_fac,1,1),Ic,SS,'un',0);

M = cell2mat_ov_sum(M_fin,xx_us,xx_uf,yy_us,yy_uf,zz_us,zz_uf,overlap_post,sizY);