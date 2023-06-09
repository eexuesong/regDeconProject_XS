% Xuesong Li 05/08/2023:

%{
Registration.m: do registration for time-sequence images of two views,
the two view images are assumed to have same pixel size and oriented 
in the same direction.

Before running the code, users should organize all the input data as TIFF 
stacks (Fig. C1): StackA_N, StackB_N … etc, where the capital letters 
A and B indicate the two perpendicular views; and the number N indicates 
the relevant time point. Here we will register B view to A view.
%}

clear all;
warning('off', 'all');

%% Load data
% load raw data
filename_data = 'StackA_0.tif';
path_data = 'D:\Code\Matlab_Code\Code_from_Min\regDeconProject-master_XS\RegistrationFusion\DataForTest';

% set parameters
t1 = 0;
t2 = 1;

%% Load DLL lib
disp('Initialize processing...');
libPath = '..\cudaLib\bin\win\';
libName = 'libapi';
libFile = strcat(libPath, libName, '.dll');
libHFile = strcat(libPath, libName, '.h');
loadlibrary(libFile, libHFile);

%% Create an output folder
path_output = strcat(path_data, '\', 'results\');
mkdir(path_output);

%% Read images
filename_A = strcat(path_data, '\', 'StackA_', num2str(t1), '.tif');
filename_B = strcat(path_data, '\', 'StackB_', num2str(t1), '.tif');
[stackA, ~] = ImageJ_formatted_TIFF.ReadTifStack(filename_A);
[stackB, ~] = ImageJ_formatted_TIFF.ReadTifStack(filename_B);
stackA = single(stackA);
stackB = single(stackB);
sizeA = size(stackA);
sizeB = size(stackB);

%% Create arguments
% results
regB = zeros(sizeA);
h_regB = libpointer('singlePtr', regB);     % registration feedback pointer: registered B stack

% input transform matrix
iTmx = eye(4);  % initial matrix for registration between two views. Tmx = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
tmxPtr = libpointer('singlePtr', iTmx);     % input matrix pointer

% input images
h_stackA = libpointer('singlePtr', stackA);     % image pointer
h_stackB = libpointer('singlePtr', stackB);     % image pointer
h_stackA_size = libpointer('uint32Ptr', sizeA); % image size pointer
h_stackB_size = libpointer('uint32Ptr', sizeB); % image size pointer

% configurations
regChoice = 2;  % *** registration choice: regChoice
                % 0: no phasor or affine registration; if flagTmx is true, transform d_img2 based on input matrix;
                % 1: phasor registraion (pixel-level translation only);
                % 2: affine registration (with or without input matrix);
                % 3: phasor registration --> affine registration (input matrix disabled);
                % 4: 2D MIP registration --> affine registration (input matrix disabled);

affMethod = 7;  % affine registration method: only if regChoice == 2, 3, 4
                % 0: no registration, but perform affine trans;
                % 1: translation only;
                % 2: rigid body; 
                % 3: 7 degrees of freedom (translation, rotation, scaling equally in 3 dimemsions)
                % 4: 9 degrees of freedom (translation, rotation, scaling); 
                % 5: 12 degrees of freedom; 
                % 6: rigid body first, then 12 degrees of freedom
                % 7: translation, rigid body, 9 degrees of freedom and then 12 degrees of freedom

flagTmx = 0;    % whether or not use input matrix, 1: yes; 0: no
                % *** flagTmx: only if regChoice == 0, 2
FTOL = 0.0001;  % reg threshold
itLimit = 2000; % maximun iteration number for registration

% parameters
deviceNum = 0;  % GPU device: numbering from 0 by CUDA;
gpuMemMode = 1; % 1: efficient GPU mode; 2: GPU memory-saved mode
verbose = 0;    % show details during processing: does not work for MATLAB
records = zeros(1, 11);
h_records = libpointer('singlePtr', records);   % reg records and feedback
tic;

%% Registration
disp('Start registration...');
for imgNum = t1:t2
    cTime1 = toc;
    disp(append('...Processing image #: ', num2str(imgNum)));

    % Read images
    filename_A = strcat(path_data, '\', 'StackA_', num2str(imgNum), '.tif');
    filename_B = strcat(path_data, '\', 'StackB_', num2str(imgNum), '.tif');
    [stackA, ~] = ImageJ_formatted_TIFF.ReadTifStack(filename_A);
    [stackB, headerB] = ImageJ_formatted_TIFF.ReadTifStack(filename_B);
    stackA = single(stackA);
    stackB = single(stackB);
    h_stackA.Value = stackA;
    h_stackB.Value = stackB;

    %         if imgNum ~= t1
    %             regChoice = 2;      % affine registration with input matrix
    %             flagTmx = 1;        % use last registration matrix as input
    %             if affMethod == 7
    %                 affMethod = 5;  % change to directly 12 DOF
    %             end
    %         end

    % run registration function: reg3d
    runStatus = calllib(libName, 'reg3d', h_regB, tmxPtr, h_stackA, h_stackB, h_stackA_size, h_stackB_size,...
            regChoice, affMethod, flagTmx, FTOL, itLimit,...
            deviceNum, gpuMemMode, verbose, h_records);
    stackB_reg = reshape(h_regB.Value, sizeA);

    % Write registered images
    if headerB.BitsPerSample == 16
        if max(stackB_reg, [], 'all') <= 65535
            stackB_reg = uint16(stackB_reg);
        else
            stackB_reg = uint16(65535 * stackB_reg ./ max(stackB_reg, [], 'all'));
        end
    end

    filename_reg = strcat(path_output, 'StackB_reg_', num2str(imgNum), '.tif');
    if isempty(headerB.resolution)
        ImageJ_formatted_TIFF.WriteTifStack(stackB_reg, filename_reg);
    elseif isempty(headerB.spacing)
        switch headerB.unit
            case 'um'
                ImageJ_formatted_TIFF.WriteTifStack(stackB_reg, filename_reg, headerB.resolution);
            case 'nm'
                ImageJ_formatted_TIFF.WriteTifStack(stackB_reg, filename_reg, headerB.resolution / 1000);
            otherwise
                ImageJ_formatted_TIFF.WriteTifStack(stackB_reg, filename_reg, headerB.resolution);
        end
    else
        switch headerB.unit
            case 'um'
                ImageJ_formatted_TIFF.WriteTifStack(stackB_reg, filename_reg, headerB.resolution, headerB.spacing);
            case 'nm'
                ImageJ_formatted_TIFF.WriteTifStack(stackB_reg, filename_reg, headerB.resolution / 1000, headerB.spacing / 1000);
            otherwise
                ImageJ_formatted_TIFF.WriteTifStack(stackB_reg, filename_reg, headerB.resolution, headerB.spacing);
        end
    end

    cTime2 = toc;
    disp(append('... ... Time cost for current image: ', num2str(cTime2 - cTime1), ' s'));
end

%% Unload DLL lib
unloadlibrary(libName);
cTime3 = toc;
disp(append('... Total time cost: ', num2str(cTime3), ' s'));
disp('Registration completed !!!');
