function I = shift_reconstruct(Y,shifts,diffphase,us_fac,Nr,Nc,Np,method,add_value)

% applies 3-d sub-pixel shifts to an input image
% INPUTS:
% Y:            input image (double 3d tensor) in space (real) or frequency (complex) domain
% shifts:       shifts 
% diffphase:    phase difference
% us_fac:       upsampling factor for subpixel shifts
% method:       method for treating boundaries
% add_value:    value to add when zero-ing out boundaries

% OUTPUT:
% I:        output image

% Written by Eftychios A. Pnevmatikakis, Simons Foundation, 2016

if isreal(Y);
    buf2ft = fftn(Y);
else
    buf2ft = Y;
end

if nargin < 9 || isempty(add_value); add_value = 0; end
if nargin < 8 || isempty(method); method = 'zero'; end
if nargin < 4 || isempty(us_fac); us_fac = 50; end

[nr,nc,np]=size(Y);
if any(shifts)
    if nargin < 5 || isempty(Nr); Nr = ifftshift(-fix(nr/2):ceil(nr/2)-1); end
    if nargin < 6 || isempty(Nc); Nc = ifftshift(-fix(nc/2):ceil(nc/2)-1); end
    if nargin < 7 || isempty(Np); Np = ifftshift(-fix(np/2):ceil(np/2)-1); end

    %shifts = shiftdim(shifts,3);
    
    row_shift = shifts(1);
    col_shift = shifts(2);
    if ismatrix(Y); pln_shift = 0; else pln_shift = shifts(3); end
    shifts = [row_shift,col_shift,pln_shift];

    if us_fac > 0
        if isvector(Nc); [Nc,Nr,Np] = meshgrid(Nc,Nr,Np); end
        Greg = buf2ft.*exp(1i*2*pi*(-row_shift*Nr/nr-col_shift*Nc/nc-pln_shift*Np/np));
    elseif us_fac == 0
        Greg = buf2ft;
    end
    Greg = Greg*exp(1i*diffphase);
    I = real(ifftn(Greg));
    I = remove_boundaries(I,shifts,method,add_value);    
else
    if isreal(Y)
        I = Y;
    else
        I = real(ifftn(Y));
    end
end

end