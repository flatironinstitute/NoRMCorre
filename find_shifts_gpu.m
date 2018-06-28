function [shifts,minY,maxY] = find_shifts_gpu(Yt,fftTempMat,fftTemp,patches,options)

fftYMat = fft(fft(fft(Yt,[],1),[],2),[],3);
Yt_patches = split_frame(Yt, patches);
fftYPatches = fft(fft(fft(Yt_patches,[],1),[],2),[],3);
shifts_gpu = register_frame_gpu(fftYMat,fftYPatches,fftTempMat,fftTemp,patches,options);
shifts_temp = gather(shifts_gpu);
shifts.shifts = shifts_temp;
shifts.shifts_up = shifts_temp;
shifts.diff = [];
minY = min(Yt(:));
maxY = max(Yt(:));