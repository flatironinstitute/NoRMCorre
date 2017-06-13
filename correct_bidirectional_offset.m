function [M,col_shift] = correct_bidirectional_offset(Y,us)

% look and correct for phase artifacts due to bi-directional scanning
% with thanks to Fabian Sinz for useful discussions.

if ~exist('us','var'); us = 10; end

sizY = size(Y);
Y = single(Y);
Y = reshape(Y,[sizY(1:2),prod(sizY(3:end))]);
mY = mean(Y,3);
mY = mY - min(mY(:));
mY1 = mY(1:2:sizY(1)-mod(sizY(1),2),:);
mY2 = mY(2:2:sizY(1),:);

mY1_us = imresize(mY1,[size(mY1,1),size(mY1,2)*us]);
mY2_us = imresize(mY2,[size(mY2,1),size(mY2,2)*us]);

fft1 = fft(mY1_us,[],2);
fft2 = fft(mY2_us,[],2);

buf_prod = fft1.*conj(fft2);

CC = ifft(buf_prod,[],2);
CCabs = abs(CC);
[~,ind] = max(mean(CCabs));
Nc_us = ifftshift(-fix(size(mY1_us,2)/2):ceil(size(mY1_us,2)/2)-1);
col_shift = Nc_us(ind)/us;

Nc = ifftshift(-fix(size(mY1,2)/2):ceil(size(mY1,2)/2)-1);

min_value = min(Y(:));
Y = Y - min_value;  % make data non-negative

Y1 = Y(1:2:end,:,:);
Y2 = Y(2:2:end,:,:);
Ys2 = real(ifft(fft(Y2,[],2).*repmat(exp(-1i*2*pi*col_shift*Nc/sizY(2)),[size(mY2,1),1,prod(sizY(3:end))]),[],2));
Ys2(:,(1:abs(round(col_shift)))*sign(col_shift) + (sizY(2)+1)*(col_shift<0),:) = Y1(1:end-mod(sizY(1),2),(1:abs(round(col_shift)))*sign(col_shift) + (sizY(2)+1)*(col_shift<0),:);
clear Y2;
M = kron(reshape(Y1,size(Y1,1),[]),[1;0]) + kron(reshape(Ys2,size(Ys2,1),[]),[0;1]);
clear Y1 Ys1;
if size(M,1) > sizY(1)
    M(sizY(1)+1:end,:) = [];
end
M = reshape(M,sizY) + min_value;