function [imgReg, oTmx, records]= reg3d_CUDA(img1, img2, libPath, libName, ...
    regChoice, affMethod, flagTmx, iTmx, FTOL, itNumLimit, devNum, gpuMemMode)
%{
GPU-based registration for 3D images
Jan 17, 2020, Min Guo

****************** Output ******************
imgReg: registered image, same size with img1
oTmx: output registration matrix
records: 11 element array, MATLAB index = C index + 1;
        [0]: actual gpu memory mode
        [1] -[3]: initial ZNCC (zero-normalized cross-correlation, negtive 
        of the cost function), intermediate ZNCC, optimized ZNCC;
        [4] -[7]: single sub iteration time (in ms), total number of sub 
        iterations, iteralation time (in s), whole registration time (in s);
        [8] -[10]: initial GPU memory, before registration, after 
        processing ( all in MB), if use gpu

****************** Input ******************
img1: target (or fixed) image
img2: source (or floating) image
libPath: folder path to the CUDA dll
libName: file name of the dll and h files 
         (the lib and h file should have same name).
regChoice: registration choice, default: 2
        0: no phasor or affine registration; if flagTmx is true, 
          transform d_img2 based on input matrix;
        1: phasor registraion (pixel-level translation only);
        2: affine registration (with or without input matrix);
        3: phasor registration --> affine registration (input matrix disabled);
        4: 2D MIP registration --> affine registration (input matrix disabled);
affMethod: affine registration method: only if regChoice == 2, 3, 4
              default: 7
        0: no registration; 
        1: translation only; 
        2: rigid body, 6 degrees of freedom (DOF); 
        3: 7 DOF (translation, rotation, scaling equally in 3 dimemsions)
        4: 9 DOF (translation, rotation, scaling); 
        5: 12 DOF; 
        6: rigid body first, then 12 DOF
        7: translation, rigid body, 9 DOF and then 12 DOF
flagTmx: whether or not use input matrix, 1: yes; 0: no
        *** flagTmx: only if regChoice == 0, 2
        true: use iTmx as input matrix;
        false: default;
iTmx: input matrix
FTOL: termination threshold, default: 0.0001
itNumLimit: maximun iteration number, default: 3000
devNum: GPU device number (CUDA convensiion, 0-based, different from
        MATLAB's 1-based, default: 0)
gpuMemMode: GPU memory mode
        1: efficient GPU mode
        2: GPU memory saved mode
%}


%% Default values;
regChoiceDefault = 2;
affMethodDefault = 7;
flagTmxDefault = 0;
iTmxDefault = eye(4);
FTOLDefault = 0.0001;
itNumLimitDefault = 3000;
devNumDefault = 0;
gpuMemModeDefault = 1;
flagVerbose = 0; % verbose not work for MATLAB
switch(nargin)
    case 4
        regChoice = regChoiceDefault;
        affMethod = affMethodDefault;
        flagTmx = flagTmxDefault;
        iTmx = iTmxDefault;
        FTOL = FTOLDefault;
        itNumLimit = itNumLimitDefault;
        devNum = devNumDefault;
        gpuMemMode = gpuMemModeDefault;
    case 5
        affMethod = affMethodDefault;
        flagTmx = flagTmxDefault;
        iTmx = iTmxDefault;
        FTOL = FTOLDefault;
        itNumLimit = itNumLimitDefault;
        devNum = devNumDefault;
        gpuMemMode = gpuMemModeDefault;
    case 6
        flagTmx = flagTmxDefault;
        iTmx = iTmxDefault;
        FTOL = FTOLDefault;
        itNumLimit = itNumLimitDefault;
        devNum = devNumDefault;
        gpuMemMode = gpuMemModeDefault;
    case 8
        FTOL = FTOLDefault;
        itNumLimit = itNumLimitDefault;
        devNum = devNumDefault;
        gpuMemMode = gpuMemModeDefault;
    case 10
        devNum = devNumDefault;
        gpuMemMode = gpuMemModeDefault;
    case 12
        % disp('all parameter manually configured')
    otherwise
        error('reg3d_CUDA: unmatched arguments, please input 4, 5, 6, 8, 10 or 12 arguments!')
end

%% load dynamic-link library
lib_load(libPath, libName);
imSize1 = size(img1);
imSize2 = size(img2);

%% Create arguments
% results 
imgReg = zeros(imSize1);
h_imgReg = libpointer('singlePtr', imgReg); % registration feedback pointer: registered image

% input images
h_img1 = libpointer('singlePtr', img1);
h_img2 = libpointer('singlePtr', img2);
    
h_imSize1 = libpointer('uint32Ptr', imSize1);   % image size pointer
h_imSize2 = libpointer('uint32Ptr', imSize2);   % image size pointer

tmxPtr = libpointer('singlePtr', iTmx); % input matrix pointer
records = zeros(1, 11);
regRecords = libpointer('singlePtr', records);  % reg records and feedback

%% Registration
cudaStatus = calllib(libName, 'reg3d', h_imgReg, tmxPtr, h_img1, h_img2, ...
    h_imSize1, h_imSize2, regChoice, affMethod, ...
    flagTmx, FTOL, itNumLimit, devNum, gpuMemMode, flagVerbose, regRecords);

if cudaStatus ~= 0
    disp('registration is probably wrong')
end

imgReg = reshape(h_imgReg.Value, imSize1);
oTmx = tmxPtr.Value;

clear h_imgReg h_img1 h_img2;

% feed back record if necessary
records = regRecords.Value;
% ncc = records(4);
% gMemPost = records(11);

% unload library
lib_unload(libName);