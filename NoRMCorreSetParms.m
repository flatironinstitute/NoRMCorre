function options = NoRMCorreSetParms(varargin)

% Struct for setting the NoRMCorre algorithm parameters. Any parameter that is
% not set gets a default value

% Author: Eftychios A. Pnevmatikakis
%            Simons Foundation, 2016

Names = [
    % dataset info
    'd1                 ' % number of rows
    'd2                 ' % number of cols
    'd3                 ' % number of planes (for 3d imaging, default: 1)
    % patches
    'grid_size          ' % size of non-overlapping regions (default: [d1,d2,d3])
    'overlap_pre        ' % size of overlapping region (default: [32,32,16])
    'min_patch_size     ' % minimum size of patch (default: [32,32,16])    
    'us_fac             ' % upsampling factor for subpixel registration (default: 20)
    'mot_uf             ' % degree of patches upsampling (default: [4,4,1])
    'max_dev            ' % maximum deviation of patch shift from rigid shift (default: [3,3,1])
    'overlap_post       ' % size of overlapping region after upsampling (default: [32,32,16])
    'max_shift          ' % maximum rigid shift in each direction (default: [15,15,5])
    'phase_flag         ' % flag for using phase correlation (default: false)
    'shifts_method      ' % method to apply shifts ('FFT','cubic','linear')
    % template updating
    'upd_template       ' % flag for online template updating (default: true)
    'init_batch         ' % length of initial batch (default: 100)
    'bin_width          ' % width of each bin (default: 10)
    'buffer_width       ' % number of local means to keep in memory (default: 50)
    'method             ' % method for averaging the template (default: {'median';'mean})
    'iter               ' % number of data passes (default: 1)
    'boundary           ' % method of boundary treatment 'NaN','copy','zero','template' (default: 'copy')
    % misc
    'add_value          ' % add dc value to data (default: 0)
    'use_parallel       ' % for each frame, update patches in parallel (default: false)
    'memmap             ' % flag for saving memory mapped motion corrected file (default: false)
    'mem_filename       ' % name for memory mapped file (default: 'motion_corrected.mat')
    'mem_batch_size     ' % batch size during memory mapping for speed (default: 5000)
    % plotting
    'plot_flag          ' % flag for plotting results in real time (default: false)
    'make_avi           ' % flag for making movie (default: false)
    'name               ' % name for movie (default: 'motion_corrected.avi')
    'fr                 ' % frame rate for movie (default: 30)
    % output type
    'output_type        ' % 'mat' (load in memory), 'memmap', 'tiff', 'hdf5', 'bin' (default:mat)
    'h5_groupname       ' % name for hdf5 dataset (default: 'mov')
    'h5_filename        ' % name for hdf5 saved file (default: 'motion_corrected.h5')
    'tiff_filename      ' % name for saved tiff stack (default: 'motion_corrected.tif')
    % use windowing
    'use_windowing      ' % flag for windowing data before fft (default: false)
    'window_length      ' % length of window on each side of the signal as a fraction of signal length
                           %    total length = length(signal)(1 + 2*window_length). (default: 0.5)
    % bitsize for reading .raw files
    'bitsize            ' % (default: 2 (uint16). other choices 1 (uint8), 4 (single), 8 (double))
    % offset from bidirectional sampling
    'correct_bidir      ' % check for offset due to bidirectional scanning (default: true)
    'nFrames            ' % number of frames to average (default: 50)
    'bidir_us           ' % upsampling factor for bidirectional sampling (default: 10)
    'col_shift          ' % known bi-directional offset provided by the user (default: [])
   ]; 
   
[m,n] = size(Names);
names = lower(Names);

% Combine all leading options structures o1, o2, ... in l1Set(o1,o2,...).
options = [];
for j = 1:m
    eval(['options.' Names(j,:) '= [];']);
end
i = 1;
while i <= nargin
    arg = varargin{i};
    if ischar(arg), break; end
    if ~isempty(arg)                      % [] is a valid options argument
        if ~isa(arg,'struct')
            error(sprintf(['Expected argument %d to be a string parameter name ' ...
                'or an options structure\ncreated with OPTIMSET.'], i));
        end
        for j = 1:m
            if any(strcmp(fieldnames(arg),deblank(Names(j,:))))
                eval(['val = arg.' Names(j,:) ';']);
            else
                val = [];
            end
            if ~isempty(val)
                eval(['options.' Names(j,:) '= val;']);
            end
        end
    end
    i = i + 1;
end

% A finite state machine to parse name-value pairs.
if rem(nargin-i+1,2) ~= 0
    error('Arguments must occur in name-value pairs.');
end
expectval = 0;                          % start expecting a name, not a value
while i <= nargin
    arg = varargin{i};
    
    if ~expectval
        if ~ischar(arg)
            error(sprintf('Expected argument %d to be a string parameter name.', i));
        end
        
        lowArg = lower(arg);
        j = strmatch(lowArg,names);
        if isempty(j)                       % if no matches
            error(sprintf('Unrecognized parameter name ''%s''.', arg));
        elseif length(j) > 1                % if more than one match
            % Check for any exact matches (in case any names are subsets of others)
            k = strmatch(lowArg,names,'exact');
            if length(k) == 1
                j = k;
            else
                msg = sprintf('Ambiguous parameter name ''%s'' ', arg);
                msg = [msg '(' deblank(Names(j(1),:))];
                for k = j(2:length(j))'
                    msg = [msg ', ' deblank(Names(k,:))];
                end
                msg = sprintf('%s).', msg);
                error(msg);
            end
        end
        expectval = 1;                      % we expect a value next
        
    else
        eval(['options.' Names(j,:) '= arg;']);
        expectval = 0;
        
    end
    i = i + 1;
end

if expectval
    error(sprintf('Expected value for parameter ''%s''.', arg));
end

Values = [
 % dataset info
    {[]} 
    {[]}
    {1}
    % patches
    {[]}                  % size of non-overlapping regions (default: [d1,d2,d3])
    {[32,32,16]}          % size of overlapping region (default: [32,32,16])
    {[32,32,16]}          % minimum size of patch (default: [32,32,16])    
    {50}                  % upsampling factor for subpixel registration (default: 50)
    {[4,4,1]}             % degree of patches upsampling (default: [4,4,1])
    {[3,3,1]}             % maximum deviation of patch shift from rigid shift (default: [3,3,1])
    {[32,32,16]}          % size of overlapping region after upsampling (default: [32,32,16])
    {[15,15,5]}           % maximum rigid shift in each direction
    {false}               % use phase correlation (good for high SNR)
    {'FFT'}               % method for applying shifts ('FFT', 'linear', 'cubic')
    % template updating
    {true}                % flag for online template updating (default: true)
    {100}                 % length of initial batch (default: 100)
    {50}                  % width of each bin (default: 10)
    {50}                  % number of local means to keep in memory (default: 50)
    {{'median';'mean'}}   % method for averaging the template (default: {'median';'mean'}
    {1}                   % number of data passes (default: 1)
    {'copy'}              % method of boundary treatment (default: 'copy')
    % misc
    {0}                   % add dc value to data (default: 0)
    {false}               % for each frame, update patches in parallel (default: false)
    {false}               % flag for saving memory mapped motion corrected file (default: false)
    {'motion_corrected.mat'} % name for memory mapped file (default: 'motion_corrected.mat')
    {1000}                % batch size used during memory mapping for faster mapping
    % plotting
    {false}               % flag for plotting results in real time (default: false)
    {false}               % flag for making movie (default: false)
    {'motion_corrected.avi'} % name for movie (default: 'motion_corrected.avi')
    {30} % frame rate for movie (default: 30)   
    % output_type    
    {'mat'}
    {'mov'}
    {'motion_corrected.h5'}
    {'motion_corrected.tif'}
    % use_windowing
    {false}
    {0.5}
    % bitsize for reading .raw files
    {2}
    % offset from bidirectional sampling
    {true}
    {50}
    {10}
    {[]}
    ];

for j = 1:m
    if eval(['isempty(options.' Names(j,:) ')'])
        eval(['options.' Names(j,:) '= Values{j};']);
    end
end

if isempty(options.d1); options.d1 = input('What is the total number of rows? \n'); end
if isempty(options.d2); options.d2 = input('What is the total number of columns? \n'); end
%if options.d3 == 1; nd = 2; else nd = 3; end
if isempty(options.grid_size); options.grid_size = [options.d1,options.d2,options.d3]; end
if length(options.grid_size) == 1; options.grid_size = options.grid_size*ones(1,3); end
if length(options.grid_size) == 2; options.grid_size(3) = 1; end
if length(options.overlap_pre) == 1; options.overlap_pre = options.overlap_pre*ones(1,3); end
if length(options.overlap_pre) == 2; options.overlap_pre(3) = 1; end
if length(options.overlap_post) == 1; options.overlap_post = options.overlap_post*ones(1,3); end
if length(options.overlap_post) == 2; options.overlap_post(3) = 1; end
if length(options.max_shift) == 1; options.max_shift = options.max_shift*ones(1,3); end
if length(options.max_shift) == 2; options.max_shift(3) = 1; end
if length(options.max_dev) == 1; options.max_dev = options.max_dev*ones(1,3); end
if length(options.max_dev) == 2; options.max_dev(3) = 1; end
if length(options.mot_uf) == 1; options.mot_uf = options.mot_uf*ones(1,3); end
if length(options.mot_uf) == 2; options.mot_uf(3) = 1; end
options.mot_uf(options.grid_size >= [options.d1,options.d2,options.d3]) = 1;