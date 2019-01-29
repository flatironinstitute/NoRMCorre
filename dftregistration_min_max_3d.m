function [output, Greg] = dftregistration_min_max_3d(buf1ft,buf2ft,usfac,min_shift,max_shift,phase_flag)
% function [output Greg] = dftregistration(buf1ft,buf2ft,usfac);
% Efficient subpixel image registration by crosscorrelation. This code
% gives the same precision as the FFT upsampled cross correlation in a
% small fraction of the computation time and with reduced memory 
% requirements. It obtains an initial estimate of the crosscorrelation peak
% by an FFT and then refines the shift estimation by upsampling the DFT
% only in a small neighborhood of that estimate by means of a 
% matrix-multiply DFT. With this procedure all the image points are used to
% compute the upsampled crosscorrelation.
% Manuel Guizar - Dec 13, 2007
%
% Rewrote all code not authored by either Manuel Guizar or Jim Fienup
% Manuel Guizar - May 13, 2016
%
% Modified by Eftychios A. Pnevmatikakis to include upper bound on possible
% shifts - November 1, 2016
%
% Citation for this algorithm:
% Manuel Guizar-Sicairos, Samuel T. Thurman, and James R. Fienup, 
% "Efficient subpixel image registration algorithms," Opt. Lett. 33, 
% 156-158 (2008).
%
% Inputs
% buf1ft    Fourier transform of reference image, 
%           DC in (1,1)   [DO NOT FFTSHIFT]
% buf2ft    Fourier transform of image to register, 
%           DC in (1,1) [DO NOT FFTSHIFT]
% usfac     Upsampling factor (integer). Images will be registered to 
%           within 1/usfac of a pixel. For example usfac = 20 means the
%           images will be registered within 1/20 of a pixel. (default = 1)
% min_shift Minimum shift in each dimension (3x1 vector). (default = -Inf, no min)
% max_shift Maximum shift in each dimension (3x1 vector). (default = Inf, no max)
%
% Outputs
% output =  [error,diffphase,net_row_shift,net_col_shift]
% error     Translation invariant normalized RMS error between f and g
% diffphase     Global phase difference between the two images (should be
%               zero if images are non-negative).
% net_row_shift net_col_shift   Pixel shifts between images
% Greg      (Optional) Fourier transform of registered version of buf2ft,
%           the global phase difference is compensated for.
%
%
% Copyright (c) 2016, Manuel Guizar Sicairos, James R. Fienup, University of Rochester
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
%     * Neither the name of the University of Rochester nor the names
%       of its contributors may be used to endorse or promote products derived
%       from this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

if ~exist('phase_flag','var')
    phase_flag = true;
end

if ~exist('usfac','var')
    usfac = 1;
end

if ~exist('max_shift','var')
    max_shift = Inf(1,3);
end

if ~exist('min_shift','var');
    min_shift = -max_shift;
end

if isscalar(max_shift); max_shift = max_shift*[1,1,1]; end
if isscalar(min_shift); min_shift = min_shift*[1,1,1]; end

[nr,nc,np]=size(buf2ft);
Nr = ifftshift(-fix(nr/2):ceil(nr/2)-1);
Nc = ifftshift(-fix(nc/2):ceil(nc/2)-1);
Np = ifftshift(-fix(np/2):ceil(np/2)-1);

buf_prod = buf1ft.*conj(buf2ft);
if usfac == 0
    % Simple computation of error and phase difference without registration
    CCmax = sum(buf1ft(:).*conj(buf2ft(:)));
    row_shift = 0;
    col_shift = 0;
    pln_shift = 0;
elseif usfac == 1
    % Single pixel registration
    if phase_flag
        buf_prod = buf_prod./abs(buf_prod);
    end
    CC = ifftn(buf_prod);
    CCabs = abs(CC);
    [~,ind] = max(CCabs(:));
    [row_shift, col_shift, pln_shift] = ind2sub([nr,nc,np],ind);
    %[row_shift, col_shift] = find(CCabs == max(CCabs(:)));
    if Nr(row_shift) > max_shift(1) || Nc(col_shift) > max_shift(2) || Np(pln_shift) > max_shift(3) || Nr(row_shift) < min_shift(1) || Nc(col_shift) < min_shift(2) || Np(pln_shift) < min_shift(3)
        CCabs2 = CCabs;
        CCabs2(Nr>max_shift(1),:,:) = 0;
        CCabs2(:,Nc>max_shift(2),:) = 0;
        CCabs2(:,:,Np>max_shift(3)) = 0;
        CCabs2(Nr<min_shift(1),:,:) = 0;
        CCabs2(:,Nc<max_shift(2),:) = 0;
        CCabs2(:,:,Np<max_shift(3)) = 0;        
        [~,ind] = max(CCabs2(:));
        [row_shift, col_shift, pln_shift] = ind2sub([nr,nc,np],ind);
        %[row_shift, col_shift] = find(CCabs == max(CCabs2(:)),1,'first');
    end    
    CCmax = CC(row_shift,col_shift,pln_shift)*nr*nc*np;
    % Now change shifts so that they represent relative shifts and not indices
    row_shift = Nr(row_shift);
    col_shift = Nc(col_shift);
    pln_shift = Np(pln_shift);
elseif usfac > 1
    % Start with usfac == 2
    buf_pad = FTpad3d(buf_prod,[2*nr,2*nc,2*np]);
    if phase_flag
        buf_pad = buf_pad./(abs(buf_pad)+1e-10);
    end
    CC = ifftn(buf_pad);
    CCabs = abs(CC);
    [~,ind] = max(CCabs(:));
    [row_shift, col_shift, pln_shift] = ind2sub([2*nr,2*nc,2*np],ind);
    %[row_shift, col_shift] = find(CCabs == max(CCabs(:)),1,'first');        
    % Now change shifts so that they represent relative shifts and not indices
    Nr2 = ifftshift(-fix(nr):ceil(nr)-1);
    Nc2 = ifftshift(-fix(nc):ceil(nc)-1);
    Np2 = ifftshift(-fix(np):ceil(np)-1);
    if Nr2(row_shift)/2 > max_shift(1) || Nc2(col_shift)/2 > max_shift(2) || Np2(pln_shift)/2 > max_shift(3) || Nr2(row_shift)/2 < min_shift(1) || Nc2(col_shift)/2 < min_shift(2) || Np2(pln_shift)/2 < min_shift(3)
        CCabs2 = CCabs;
        CCabs2(Nr2/2>max_shift(1),:,:) = 0;
        CCabs2(:,Nc2/2>max_shift(2),:) = 0;
        CCabs2(:,:,Np2/2>max_shift(3)) = 0;
        CCabs2(Nr2/2<min_shift(1),:,:) = 0;
        CCabs2(:,Nc2/2<min_shift(2),:) = 0;
        CCabs2(:,:,Np2/2<min_shift(3)) = 0;
        [~,ind] = max(CCabs2(:));
        [row_shift, col_shift, pln_shift] = ind2sub([2*nr,2*nc,2*np],ind);
        %CCabs2 = CCabs(abs(Nr)<=max_shift(1),abs(Nc)<=max_shift(2));
        %[row_shift, col_shift] = find(CCabs == max(CCabs2(:)),1,'first');
    end     
    CCmax = CC(row_shift,col_shift,pln_shift)*nr*nc*np;
    row_shift = Nr2(row_shift)/2;
    col_shift = Nc2(col_shift)/2;
    pln_shift = Np2(pln_shift)/2;
    % If upsampling > 2, then refine estimate with matrix multiply DFT
    if usfac > 2,
        %%% DFT computation %%%
        % Initial shift estimate in upsampled grid
        row_shift = round(row_shift*usfac)/usfac; 
        col_shift = round(col_shift*usfac)/usfac;
        pln_shift = round(pln_shift*usfac)/usfac;
        dftshift = fix(ceil(usfac*1.5)/2); %% Center of output array at dftshift+1
        % Matrix multiply DFT around the current shift estimate
        CC = conj(dftups3d(buf2ft.*conj(buf1ft),ceil(usfac*1.5),ceil(usfac*1.5),ceil(usfac*1.5),usfac,...
            dftshift-row_shift*usfac,dftshift-col_shift*usfac,dftshift-pln_shift*usfac));
        % Locate maximum and map back to original pixel grid 
        CCabs = abs(CC);
        [~,ind] = max(CCabs(:));
        [rloc, cloc, ploc] = ind2sub(size(CC),ind);
        %[rloc, cloc] = find(CCabs == max(CCabs(:)),1,'first');
        CCmax = CC(rloc,cloc,ploc);
        rloc = rloc - dftshift - 1;
        cloc = cloc - dftshift - 1;
        ploc = ploc - dftshift - 1;
        row_shift = row_shift + rloc/usfac;
        col_shift = col_shift + cloc/usfac;
        pln_shift = pln_shift + ploc/usfac;
    end

    % If its only one row or column the shift along that dimension has no
    % effect. Set to zero.
    if nr == 1
        row_shift = 0;
    end
    if nc == 1
        col_shift = 0;
    end
    if np == 1
        pln_shift = 0;
    end    
end  

if 0
    rg00 = sum(abs(buf1ft(:)).^2);
    rf00 = sum(abs(buf2ft(:)).^2);
    error = 1.0 - abs(CCmax).^2/(rg00*rf00);
    error = sqrt(abs(error));
else
    error = 0;
end
diffphase = angle(CCmax);

output=[error,diffphase,row_shift,col_shift,pln_shift];

% Compute registered version of buf2ft
if (nargout > 1)&&(usfac > 0)
    [Nc,Nr,Np] = meshgrid(Nc,Nr,Np);
    Greg = buf2ft.*exp(1i*2*pi*(-row_shift*Nr/nr-col_shift*Nc/nc-pln_shift*Np/np));
    Greg = Greg*exp(1i*diffphase);
elseif (nargout > 1)&&(usfac == 0)
    Greg = buf2ft*exp(1i*diffphase);
end
return

function out=dftups3d(in,nor,noc,nop,usfac,roff,coff,poff)
% function out=dftups(in,nor,noc,usfac,roff,coff);
% Upsampled DFT by matrix multiplies, can compute an upsampled DFT in just
% a small region.
% usfac         Upsampling factor (default usfac = 1)
% [nor,noc]     Number of pixels in the output upsampled DFT, in
%               units of upsampled pixels (default = size(in))
% roff, coff    Row and column offsets, allow to shift the output array to
%               a region of interest on the DFT (default = 0)
% Recieves DC in upper left corner, image center must be in (1,1) 
% Manuel Guizar - Dec 13, 2007
% Modified from dftus, by J.R. Fienup 7/31/06

% This code is intended to provide the same result as if the following
% operations were performed
%   - Embed the array "in" in an array that is usfac times larger in each
%     dimension. ifftshift to bring the center of the image to (1,1).
%   - Take the FFT of the larger array
%   - Extract an [nor, noc] region of the result. Starting with the 
%     [roff+1 coff+1] element.

% It achieves this result by computing the DFT in the output array without
% the need to zeropad. Much faster and memory efficient than the
% zero-padded FFT approach if [nor noc] are much smaller than [nr*usfac nc*usfac]

[nr,nc,np]=size(in);
% Set defaults
if exist('roff', 'var')~=1, roff=0;  end
if exist('coff', 'var')~=1, coff=0;  end
if exist('poff', 'var')~=1, poff=0;  end
if exist('usfac','var')~=1, usfac=1; end
if exist('noc',  'var')~=1, noc=nc;  end
if exist('nor',  'var')~=1, nor=nr;  end
if exist('nop',  'var')~=1, nop=np;  end

%tp3d =  @(x,y,z) reshape(kron(y(:),x(:))*z(:)',length(x),length(y),length(z));

% Compute kernels and obtain DFT by matrix products
%kernc=exp((-1i*2*pi/(nc*usfac))*( ifftshift(0:nc-1).' - floor(nc/2) )*( (0:noc-1) - coff ));
kernc=exp((-1i*2*pi/(nc*usfac))*((0:noc-1).' - coff) * (ifftshift(0:nc-1) - floor(nc/2)) );
kernr=exp((-1i*2*pi/(nr*usfac))*((0:nor-1).' - roff) * (ifftshift(0:nr-1) - floor(nr/2)) );
kernp=exp((-1i*2*pi/(np*usfac))*((0:nop-1).' - poff) * (ifftshift(0:np-1) - floor(np/2)) );

%kernc = exp((-1i*2*pi/(nc*usfac))*tp3d(ifftshift(0:nc-1) - floor(nc/2),(0:noc-1)-coff,(0:noc-1)-coff));
%kernr = exp((-1i*2*pi/(nr*usfac))*tp3d((0:nor-1)-roff,ifftshift(0:nr-1) - floor(nr/2),(0:nor-1)-roff));
%kernp = exp((-1i*2*pi/(np*usfac))*tp3d((0:nop-1)-poff,(0:nop-1)-poff,ifftshift(0:np-1) - floor(np/2)));
out = reshape(kernr*reshape(in,nr,[]),nor,nc,np);
out = permute(out,[2,1,3]);
out = reshape(kernc*reshape(out,nc,[]),noc,nor,np);
out = permute(out,[2,1,3]);
out = permute(out,[3,2,1]);
out = reshape(kernp*reshape(out,np,[]),noc,nor,nop);
out = permute(out,[3,2,1]);

%out=kernr*in*kernc;
return


function [ imFTout ] = FTpad3d(imFT,outsize)
% imFTout = FTpad3d(imFT,outsize)
% Pads or crops the Fourier transform to the desired ouput size. Taking 
% care that the zero frequency is put in the correct place for the output
% for subsequent FT or IFT. Can be used for Fourier transform based
% interpolation, i.e. dirichlet kernel interpolation. 
%
%   Inputs
% imFT      - Input complex array with DC in [1,1,1]
% outsize   - Output size of array [ny nx np] 
%
%   Outputs
% imout   - Output complex image with DC in [1,1,1]
% Manuel Guizar - 2014.06.02
% modified by Eftychios A. Pnevmatikakis to include 3d analysis and max
% shifts - 2016.11.03

Nout = outsize;
Nin = size(imFT);
imFT = fftshift(imFT);
center = floor(size(imFT)/2)+1;
imFTout = zeros(outsize);
centerout = floor(size(imFTout)/2)+1;

% imout(centerout(1)+[1:Nin(1)]-center(1),centerout(2)+[1:Nin(2)]-center(2)) ...
%     = imFT;
cenout_cen = centerout - center;

imFTout(max(cenout_cen(1)+1,1):min(cenout_cen(1)+Nin(1),Nout(1)),max(cenout_cen(2)+1,1):min(cenout_cen(2)+Nin(2),Nout(2)),max(cenout_cen(3)+1,1):min(cenout_cen(3)+Nin(3),Nout(3))) ...
    = imFT(max(-cenout_cen(1)+1,1):min(-cenout_cen(1)+Nout(1),Nin(1)),max(-cenout_cen(2)+1,1):min(-cenout_cen(2)+Nout(2),Nin(2)),max(-cenout_cen(3)+1,1):min(-cenout_cen(3)+Nout(3),Nin(3)));

%imFTout2 = padarray(imFT,cenout_cen,0,'pre');
%imFTout = padarray(imFTout2,cenout_cen-mod(size(imFT),2),0,'post');
imFTout = ifftshift(imFTout)*Nout(1)*Nout(2)*Nout(3)/(Nin(1)*Nin(2)*Nin(3));
return