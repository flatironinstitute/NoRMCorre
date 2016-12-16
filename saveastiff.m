function res = saveastiff(data, path, options)
% options.color
%   : true or FALSE
%   : If this is true, third dimension should be 3 and the data is saved as a color image.
% options.compress
%   : 'no', 'lzw', 'jpeg' or 'adobe'.
%     Compression type.
%       'no'    : Uncompressed(Default)
%       'lzw'   : lossless LZW
%       'jpeg'  : lossy JPEG (When using JPEG compression, ImageWidth,
%                 ImageLength, and RowsPerStrip must be multiples of 16.)
%       'adobe' : lossless Adobe-style
% options.message
%   : TRUE or false.
%     If this is false, all messages are skipped. 
% options.append
%   : true or FALSE
%     If path is exist, the data is appended to an existing file.
%     If path is not exist, this options is ignored.
% options.overwrite
%   : true or FALSE
%     Overwrite to an existing file.
% options.big 
%   : true or FALSE, 
%     Use 64 bit addressing and allows for files > 4GB
% 
% Defalut value of 'options' is
%     options.color     = false;
%     options.compress  = 'no';
%     options.message   = true;
%     options.append    = false;
%     options.overwrite = false;
%     options.big       = false;
% 
% res : Return value. It is 0 when the function is finished with no error.
%       If an error is occured in the function, it will have a positive
%       number (error code).
%
% Copyright (c) 2012, YoonOh Tak
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
%     * Neither the name of the Gwangju Institute of Science and Technology (GIST), Republic of Korea nor the names 
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

tStart = tic;
errcode = 0;
try
%% Init options parameter    
if nargin < 3 % Use default options
    options.color = false;
    options.compress = 'no';
    options.message = true;
    options.append = false;
    options.overwrite = false;
end
if ~isfield(options, 'message'),   options.message   = true; end
if ~isfield(options, 'append'),    options.append    = false; end
if ~isfield(options, 'compress'),  options.compress  = 'no';  end
if ~isfield(options, 'color'),     options.color     = false; end
if ~isfield(options, 'overwrite'), options.overwrite = false; end
if  isfield(options, 'big') == 0,  options.big       = false; end

if isempty(data), errcode = 1; assert(false); end
if (options.color == false && ndims(data) > 3) || ...
   (options.color == true && ndims(data) > 4)
    % Maximum dimension of a grayscale image is 3 of [height, width, frame]
    % Maximum dimension of a color image is 4 of [height, width, color, frame]
    errcode = 2; assert(false);
end

%% Get image informations
% http://www.awaresystems.be/imaging/tiff/tifftags/photometricinterpretation.html
if ~options.color
    if ndims(data) >= 4, errcode = 2; assert(false); end;
    [height, width, depth] = size(data);
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
%     tagstruct.Photometric = Tiff.Photometric.MinIsWhite;
%     tagstruct.Photometric = Tiff.Photometric.Mask;
%     tagstruct.Photometric = Tiff.Photometric.Separated;
else
    if ndims(data) >= 5, errcode = 2; assert(false); end;
    [height, width, cc, depth] = size(data); % cc: color channels. 3: rgb, 4: rgb with alpha channel
    if cc ~= 3 && cc ~= 4, errcode = 3; assert(false); end;
    tagstruct.Photometric = Tiff.Photometric.RGB;
%     tagstruct.Photometric = Tiff.Photometric.CIELab;
%     tagstruct.Photometric = Tiff.Photometric.ICCLab;
%     tagstruct.Photometric = Tiff.Photometric.ITULab;
%     (Unsupported)tagstruct.Photometric = Tiff.Photometric.Palette;
%     (Unsupported)tagstruct.Photometric = Tiff.Photometric.YCbCr;
end
tagstruct.ImageLength = height;
tagstruct.ImageWidth = width;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky; % (RGB RGB,RGB RGB,RGB RGB), http://www.awaresystems.be/imaging/tiff/tifftags/planarconfiguration.html
% tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Separate; % (RRR RRR, GGG GGG, BBB BBB), % http://www.awaresystems.be/imaging/tiff/tifftags/planarconfiguration.html


%% Complex number
% http://www.awaresystems.be/imaging/tiff/tifftags/samplesperpixel.html
if ~options.color && isreal(data) % Grayscale image with real numbers
    tagstruct.SamplesPerPixel = 1;
    data = reshape(data, height, width, 1, depth);
elseif ~options.color && ~isreal(data) % Grayscale image with complex numbers
    tagstruct.SamplesPerPixel = 2;
    data = reshape([real(data) imag(data)], height, width, 2, depth);
elseif options.color && isreal(data) % Color image with real numbers
    tagstruct.SamplesPerPixel = cc;
    if cc == 4
        tagstruct.ExtraSamples = Tiff.ExtraSamples.AssociatedAlpha; % The forth channel is alpha channel
    end
    data = reshape(data, height, width, cc, depth);
elseif options.color && ~isreal(data) % Color image with complex numbers
    tagstruct.SamplesPerPixel = cc * 2;
    if cc == 3
        tagstruct.ExtraSamples = repmat(Tiff.ExtraSamples.Unspecified, 1, 3); % 3(real)+3(imag) = 6 = 3(rgb) + 3(Extra)
    else
        tagstruct.ExtraSamples = repmat(Tiff.ExtraSamples.Unspecified, 1, 5); % 4(real)+4(imag) = 8 = 3(rgb) + 5(Extra)
    end
    data = reshape([real(data) imag(data)], height, width, cc*2, depth);
end

%% Image compression
% http://www.awaresystems.be/imaging/tiff/tifftags/compression.html
switch lower(options.compress)
    case 'no'
        tagstruct.Compression = Tiff.Compression.None;
    case 'lzw'
        tagstruct.Compression = Tiff.Compression.LZW;
    case 'jpeg'
        tagstruct.Compression = Tiff.Compression.JPEG;
    case 'adobe'
        tagstruct.Compression = Tiff.Compression.AdobeDeflate;
    otherwise
        % Use tag nubmer in http://www.awaresystems.be/imaging/tiff/tifftags/compression.html
        tagstruct.Compression = options.compress;
end

%% Sample format
% http://www.awaresystems.be/imaging/tiff/tifftags/sampleformat.html
switch class(data)
    % Unsupported Matlab data type: char, logical, cell, struct, function_handle, class.
    case {'uint8', 'uint16', 'uint32'}
        tagstruct.SampleFormat = Tiff.SampleFormat.UInt;
    case {'int8', 'int16', 'int32'}
        tagstruct.SampleFormat = Tiff.SampleFormat.Int;
        if options.color
            errcode = 4; assert(false);
        end
    case {'single', 'double', 'uint64', 'int64'}
        tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
    otherwise
        % (Unsupported)Void, ComplexInt, ComplexIEEEFP
        errcode = 5; assert(false);
end

%% Bits per sample
% http://www.awaresystems.be/imaging/tiff/tifftags/bitspersample.html
switch class(data)
    case {'uint8', 'int8'}
        tagstruct.BitsPerSample = 8;
    case {'uint16', 'int16'}
        tagstruct.BitsPerSample = 16;
    case {'uint32', 'int32'}
        tagstruct.BitsPerSample = 32;
    case {'single'}
        tagstruct.BitsPerSample = 32;
    case {'double', 'uint64', 'int64'}
        tagstruct.BitsPerSample = 64;
    otherwise
        errcode = 5; assert(false);
end

%% Rows per strip
maxstripsize = 8*1024;
tagstruct.RowsPerStrip = ceil(maxstripsize/(width*(tagstruct.BitsPerSample/8)*size(data,3))); % http://www.awaresystems.be/imaging/tiff/tifftags/rowsperstrip.html
if tagstruct.Compression == Tiff.Compression.JPEG
    tagstruct.RowsPerStrip = max(16,round(tagstruct.RowsPerStrip/16)*16);
end

%% Overwrite check
if exist(path, 'file') && ~options.append
    if ~options.overwrite
        errcode = 6; assert(false);
    end
end

%% Save path configuration
path_parent = pwd;
[pathstr, fname, fext] = fileparts(path);
if ~isempty(pathstr)
    if ~exist(pathstr, 'dir')
        mkdir(pathstr);
    end
    cd(pathstr);
end

%% Write image data to a file
file_opening_error_count = 0;
while ~exist('tfile', 'var')
    try
        if ~options.append % Make a new file
            s=whos('data');
            if s.bytes > 2^32-1 || options.big
                tfile = Tiff([fname, fext], 'w8'); % Big Tiff file
            else
                tfile = Tiff([fname, fext], 'w');
            end
        else
            if ~exist([fname, fext], 'file') % Make a new file
                s=whos('data');
                if s.bytes > 2^32-1 || options.big
                    tfile = Tiff([fname, fext], 'w8'); % Big Tiff file
                else
                    tfile = Tiff([fname, fext], 'w');
                end
            else % Append to an existing file
                tfile = Tiff([fname, fext], 'r+');
                while ~tfile.lastDirectory(); % Append a new image to the last directory of an exiting file
                    tfile.nextDirectory();
                end
                tfile.writeDirectory();
            end
        end
    catch
        file_opening_error_count = file_opening_error_count + 1;
        pause(0.1);
        if file_opening_error_count > 5 % automatically retry to open for 5 times.
            reply = input('Failed to open the file. Do you wish to retry? Y/n: ', 's');
            if isempty(reply) || any(upper(reply) == 'Y')
                file_opening_error_count = 0;
            else
                errcode = 7;
                assert(false);
            end
        end
    end
end

for d = 1:depth
    tfile.setTag(tagstruct);
    tfile.write(data(:, :, :, d));
    if d ~= depth
       tfile.writeDirectory();
    end
end

tfile.close();
if exist('path_parent', 'var'), cd(path_parent); end

tElapsed = toc(tStart);
if options.message
    display(sprintf('The file was saved successfully. Elapsed time : %.3f s.', tElapsed));
end

catch exception
%% Exception management
    if exist('tfile', 'var'), tfile.close(); end
    switch errcode
        case 1
            if options.message, error '''data'' is empty.'; end;
        case 2
            if options.message, error 'Data dimension is too large.'; end;
        case 3
            if options.message, error 'Third dimesion (color depth) should be 3 or 4.'; end;
        case 4
            if options.message, error 'Color image cannot have int8, int16 or int32 format.'; end;
        case 5
            if options.message, error 'Unsupported Matlab data type. (char, logical, cell, struct, function_handle, class)'; end;
        case 6
            if options.message, error 'File already exists.'; end;
        case 7
            if options.message, error(['Failed to open the file ''' path '''.']); end;
        otherwise
            if exist('fname', 'var') && exist('fext', 'var')
                delete([fname fext]);
            end
            if exist('path_parent', 'var'), cd(path_parent); end
            rethrow(exception);
    end
    if exist('path_parent', 'var'), cd(path_parent); end
end
res = errcode;
end