% Xuesong Li 05/09/2023:

%{
Fusion.m: do registration and re for time-sequence images of two views,
the two view images (input) are assumed to have isotripic pixel size and oriented 
in the same direction.

Before running the code, users should organize all the input data as TIFF 
stacks (Fig. C1): StackA_N, StackB_N … etc, where the capital letters 
A and B indicate the two perpendicular views; and the number N indicates 
the relevant time point. Here we will register B view to A view.

Users need to generate forward projectors and Wiener-Butterworth backward 
projectors according to the Online Methods and Supplementary Note 2 
(or the BackProjector.m script referenced above). The forward and backward 
projectors should be oriented from the same perspective as their 
corresponding views. These projector images can be placed in the same 
folder as the input data or grouped into another folder, named PSFA.tif 
for the forward projector of the first view, PSFA_BP.tif for the backward 
projector of the first view, PSFB.tif for the forward projector of the 
second view, and PSFB_BP.tif for the forward projector of the second view.

Note that if the backward projectors are the transpose of the forward 
projectors, then the deconvolution is equivalent to traditional 
Richardson-Lucy deconvolution.

While running the Fusion.m script, three dialog boxes will pop up to guide 
users to 
(1) upload the raw dual-view data (see Fig. C1); 
(2) upload the forward and backward projectors; 
(3) set parameters (similar to Fig. C4 except there is an input for setting
 the iteration number) to
a) choose registration mode:
    1: translation only;
    2: rigid body;
    3: 7 degrees of freedom (translation, rotation, scaling equally in 3 dimensions);
    4: 9 degrees of freedom (translation, rotation, scaling);
    5: 12 degrees of freedom (translation, rotation, scaling, shearing);
b) set the iteration number (the typical number of iterations we used for 
Wiener-Butterworth deconvolution and traditional deconvolution is 1 and 10-20, respectively);
c) set the number of time points (in range form, e.g. 0-9) to be processed.

After processing, the registered B view data and joint deconvolution 
outcomes will be saved in the same folder that contains the input data. 
They are saved as StackB_reg_N and Decon_N, where N indicates the number 
of time points.
%}

clear all;
% warning('off', 'all');

%% GUI part
% load raw data
filename_data = 'StackA_0.tif';
path_data = 'D:\Code\Matlab_Code\Code_from_Min\regDeconProject-master_XS\RegistrationFusion\DataForTest';
% load PSF data
path_psf = 'D:\Code\Matlab_Code\Code_from_Min\regDeconProject-master_XS\RegistrationFusion\DataForTest';

% set parameters
itNum = 10;
t1 = 0;
t2 = 1;

%% Load DLL lib
disp('Initialize processing...');
libPath = '..\cudaLib\bin\win\';
libName = 'libapi';

%% Create an output folder
path_output = [path_data, '\', 'results\'];
mkdir(path_output);

%% Forward and back projectors
[PSFA, ~] = ImageJ_formatted_TIFF.ReadTifStack(strcat(path_psf, '\', 'PSFA.tif'));
[PSFB, ~] = ImageJ_formatted_TIFF.ReadTifStack(strcat(path_psf, '\', 'PSFB.tif'));
[PSFA_bp, ~] = ImageJ_formatted_TIFF.ReadTifStack(strcat(path_psf, '\', 'PSFA_BP.tif'));
[PSFB_bp, ~] = ImageJ_formatted_TIFF.ReadTifStack(strcat(path_psf, '\', 'PSFB_BP.tif'));
PSFA = single(PSFA);
PSFB = single(PSFB);
PSFA_bp = single(PSFA_bp);
PSFB_bp = single(PSFB_bp);

%% Create arguments
% input transform matrix
iTmx = eye(4);  % initial matrix for registration between two views. Tmx = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];

%  configurations
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
tic;

disp('Start processing...');
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

    if imgNum ~= t1
        regChoice = 2;      % affine registration with input matrix
        flagTmx = 1;        % use last registration matrix as input
        iTmx = oTmx;
        if affMethod == 7
            affMethod = 5;  % change to directly 12 DOF
        end
    end

    %% Registration
    disp('... ... Performing registration ...');

    % run registration function: reg3d
    [stackB_reg, oTmx] = reg3d_CUDA(stackA, stackB, libPath, libName, regChoice, affMethod, flagTmx, iTmx, FTOL, itLimit, deviceNum, gpuMemMode);

    %% Deconvolution
    disp('... ... Performing deconvolution ...');
    flagConstInitial = 0;
    flagUnmatch = 0;

    % run registration function: decon_dualview
    stack_decon = decon_dualview_CUDA(stackA, stackB_reg, PSFA, PSFB, PSFA_bp, PSFB_bp, libPath, libName,...
        flagConstInitial, itNum, flagUnmatch, deviceNum, gpuMemMode);

    %% Write registered images
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

    %% Write deconvolved images
    if headerB.BitsPerSample == 16
        if max(stack_decon, [], 'all') <= 65535
            stack_decon = uint16(stack_decon);
        else
            stack_decon = uint16(65535 * stack_decon ./ max(stack_decon, [], 'all'));
        end
    end

    filename_decon = strcat(path_output, 'Decon_', num2str(imgNum), '.tif');
    if isempty(headerB.resolution)
        ImageJ_formatted_TIFF.WriteTifStack(stack_decon, filename_decon);
    elseif isempty(headerB.spacing)
        switch headerB.unit
            case 'um'
                ImageJ_formatted_TIFF.WriteTifStack(stack_decon, filename_reg, headerB.resolution);
            case 'nm'
                ImageJ_formatted_TIFF.WriteTifStack(stack_decon, filename_reg, headerB.resolution / 1000);
            otherwise
                ImageJ_formatted_TIFF.WriteTifStack(stack_decon, filename_reg, headerB.resolution);
        end
    else
        switch headerB.unit
            case 'um'
                ImageJ_formatted_TIFF.WriteTifStack(stack_decon, filename_reg, headerB.resolution, headerB.spacing);
            case 'nm'
                ImageJ_formatted_TIFF.WriteTifStack(stack_decon, filename_reg, headerB.resolution / 1000, headerB.spacing / 1000);
            otherwise
                ImageJ_formatted_TIFF.WriteTifStack(stack_decon, filename_reg, headerB.resolution, headerB.spacing);
        end
    end

    cTime2 = toc;
    disp(append('... ... Time cost for current image: ', num2str(cTime2 - cTime1), ' s'));
end

cTime3 = toc;
disp(append('... Total time cost: ', num2str(cTime3), ' s'));
disp('Processing completed !!!');

