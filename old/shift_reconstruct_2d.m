function I = shift_reconstruct_2d(Y,shifts,us_fac,zero_rows,zero_cols,Nr,Nc)

% applies 2-d sub-pixel shifts to an input image
% INPUTS:
% Y:            input image (double matrix)
% shifts:       shifts 
% us_fac:       upsampling factor for subpixel shifts
% zero_rows:    zero out shifted rows part 
% zero_cols:    zero out shifted columns part
% Nr:           precomputed meshgrid
% Nc:           precomputed meshgrid

% OUTPUT:
% I:        output image

% Written by Eftychios A. Pnevmatikakis, Simons Foundation, 2016

[nr,nc]=size(Y);
if any(shifts)
    if any(shifts - round(shifts))
        if nargin < 5 || isempty(Nr);  Nr = ifftshift(-fix(nr/2):ceil(nr/2)-1); end
        if nargin < 6 || isempty(Nc);  Nc = ifftshift(-fix(nc/2):ceil(nc/2)-1); end
        row_shift = shifts(1);
        col_shift = shifts(2);
        buf2ft = fft2(Y);

        if us_fac > 0
            if ~ismatrix(Nc) || ~ismatrix(Nr)
                [Nc,Nr] = meshgrid(Nc,Nr);
            end
            Greg = buf2ft.*exp(1i*2*pi*(-row_shift*Nr/nr-col_shift*Nc/nc));
        elseif us_fac == 0
            Greg = buf2ft;
        end

        I = abs(ifft2(Greg));
        shifts = round(shifts);
    else
        I = circshift(Y,squeeze(shifts));
    end
    if zero_rows
        I((1:abs(shifts(1)))*sign(shifts(1)) + (nr+1)*(shifts(1)<0),:) = 0;
    end

    if zero_cols
        I(:,(1:abs(shifts(2)))*sign(shifts(2)) + (nc+1)*(shifts(2)<0)) = 0;
    end
else
    I = Y;
end