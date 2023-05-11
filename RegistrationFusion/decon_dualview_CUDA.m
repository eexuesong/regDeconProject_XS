function [stack_decon, records]= decon_dualview_CUDA(stackA, stackB, PSFA, PSFB, PSFA_bp, PSFB_bp, libPath, libName, ...
    flagConstInitial, itNum, flagUnmatch, devNum, gpuMemMode)
%{
GPU-based dual-view joint deconvolution for diSPIM
Jan 17, 2020, Min Guo

****************** Output ******************
stack_decon: deconvolved image, same size with stackA
records: 11 element array, MATLAB index = C index + 1;
        [0]: actual gpu memory mode
        [1] -[3]: initial ZNCC (zero-normalized cross-correlation, negtive 
        of the cost function), intermediate ZNCC, optimized ZNCC;
        [4] -[7]: single sub iteration time (in ms), total number of sub 
        iterations, iteralation time (in s), whole registration time (in s);
        [8] -[10]: initial GPU memory, before registration, after 
        processing ( all in MB), if use gpu

****************** Input ******************
stackA: view A image
stackB: view B image
PSFA: view A forward projector
PSFB: view B forward projector
PSFA_bp: view A back projector
PSFB_bp: view B back projector
libPath: folder path to the CUDA dll
libName: file name of the dll and h files 
         (the lib and h file should have same name).
flagConstInitial: whether or not to use a const as the initial estimate, 
        default: 0
itNum: iteration number of RL deconvolution, default: 10
flagUnmatch: whether or not to use the unmatched Wiener-Butterworth back 
        projectors, default: 0
devNum: GPU device number (CUDA convensiion, 0-based, different from
        MATLAB's 1-based, default: 0)
gpuMemMode: GPU memory mode
        1: efficient GPU mode
        2: GPU memory saved mode
%}


%% Default values;
flagConstInitialDefault = 0;
itNumDefault = 10;
flagUnmatchDefault = 0;
devNumDefault = 0;
gpuMemModeDefault = 1;
flagVerbose = 0; % verbose not work for MATLAB
switch(nargin)
    case 8
        flagConstInitial = flagConstInitialDefault;
        itNum = itNumDefault;
        flagUnmatch = flagUnmatchDefault;
        devNum = devNumDefault;
        gpuMemMode = gpuMemModeDefault;
    case 9
        itNum = itNumDefault;
        flagUnmatch = flagUnmatchDefault;
        devNum = devNumDefault;
        gpuMemMode = gpuMemModeDefault;
    case 10
        flagUnmatch = flagUnmatchDefault;
        devNum = devNumDefault;
        gpuMemMode = gpuMemModeDefault;
    case 11
        devNum = devNumDefault;
        gpuMemMode = gpuMemModeDefault;
    case 13
        % disp('all parameter manually configured')
    otherwise
        error('decon_dualview_CUDA: unmatched arguments, please input 8, 9, 10, 11 or 13 arguments!')
end

%% load dynamic-link library
lib_load(libPath, libName);
stack_size = size(stackA);
PSF_size = size(PSFA);

%% Create arguments
% results
stack_decon = zeros(stack_size, 'single');
h_decon = libpointer('singlePtr', stack_decon); % decon feedback pointer: decon result

% input images
h_stackA = libpointer('singlePtr', stackA);         % image pointer
h_stackB = libpointer('singlePtr', stackB);         % image pointer
h_stack_size = libpointer('uint32Ptr', stack_size); % image size pointer

% input PSFs
h_PSFA = libpointer('singlePtr', PSFA);             % PSF pointer
h_PSFB = libpointer('singlePtr', PSFB);             % PSF pointer
h_PSFA_bp = libpointer('singlePtr', PSFA_bp);       % PSF pointer
h_PSFB_bp = libpointer('singlePtr', PSFB_bp);       % PSF pointer
h_PSF_size = libpointer('uint32Ptr', PSF_size);     % PSF size pointer

% parameters
records = zeros(1, 11);
h_records = libpointer('singlePtr', records);  % reg records and feedback

%% Registration
cudaStatus = calllib(libName, 'decon_dualview', h_decon, h_stackA, h_stackB, h_stack_size, h_PSFA, h_PSFB, h_PSF_size,...
            flagConstInitial, itNum, devNum, gpuMemMode, flagVerbose, h_records, flagUnmatch, h_PSFA_bp, h_PSFB_bp);

if cudaStatus ~= 0
    disp('deconvolution is probably wrong')
end

stack_decon = reshape(h_decon.Value, stack_size);

clear h_decon h_stackA h_stackB h_PSFA h_PSFB h_PSFA_bp h_PSFB_bp;

% feed back record if necessary
records = h_records.Value;
% ncc = records(4);
% gMemPost = records(11);

% unload library
lib_unload(libName);