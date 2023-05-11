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
[filename_data, path_data] = uigetfile('*.tif', 'Choose any one of raw data');
% load PSF data
[filename_psf, path_psf] = uigetfile('*.tif', 'Choose any one of PSF image', path_data);

% set parameters
inputDialog = uifigure('Name', 'Set Parameters', 'Position', [300 500 400 320]);
% set registration mode
uilabel(inputDialog, 'Text', 'Registration mode:', 'Position', [30 290 300 20]);
% Create a button group and radio buttons:
bg = uibuttongroup('Parent', inputDialog, 'Position', [50 170 300 120]);
rb1 = uiradiobutton(bg, 'Text', 'translation only', 'Position', [10 90 200 20]);
rb2 = uiradiobutton(bg, 'Text', 'rigid body', 'Position', [10 70 200 20]);
rb3 = uiradiobutton(bg, 'Text', '7 DOF', 'Position', [10 50 200 20]);
rb4 = uiradiobutton(bg, 'Text', '9 DOF', 'Position', [10 30 200 20]);
rb5 = uiradiobutton(bg, 'Text', '12 DOF', 'Position', [10 10 200 20]);
% for time points
uilabel(inputDialog, 'Text', 'Enter iteration number: :', 'Position', [30 140 300 20]);
txa1 = uitextarea(inputDialog, 'Value', '10', 'Position', [50 120 100 20]);
uilabel(inputDialog, 'Text', 'Enter time points to be processed:', 'Position', [30 100 300 20]);
txa2 = uitextarea(inputDialog, 'Value', '0-1', 'Position', [50 70 100 30]);
% buttons to confirm/cancel processing
pb1 = uibutton(inputDialog, 'push', 'Position', [130 30 50 20], 'Text', 'Yes',...
    'ButtonPushedFcn', @(pb1,event) pbYes(inputDialog, bg, txa1, txa2, path_data, path_psf));
pb2 = uibutton(inputDialog, 'push', 'Position', [210 30 50 20], 'Text', 'Cancel',...
    'ButtonPushedFcn', @(pb1,event) pbCancel(inputDialog));

%% Cancel dialog box
function pbCancel(inputDialog)
    delete(inputDialog);
    disp('Processing cancelled!!!')
end

%% Yes dialog box
function pbYes(inputDialog, bg, txa1, txa2, path_data, path_psf)
    bgChoice = get(get(bg, 'SelectedObject'), 'Text');
    itStr = get(txa1, 'Value');
    itNum = str2double(itStr{1});
    timeStr = get(txa2,'Value');
    timepoints = strsplit(timeStr{1}, '-');
    t1 = str2double(timepoints{1});
    t2 = str2double(timepoints{2});
    delete(inputDialog);
    disp('Initialize processing...');

    %% Load DLL lib
    libPath = '..\cudaLib\bin\win\';
    libName = 'libapi';

    %% Create an output folder
    path_output = [path_data, 'results\'];
    mkdir(path_output);
    
    %% Forward and back projectors
    [PSFA, ~] = ImageJ_formatted_TIFF.ReadTifStack(strcat(path_psf, 'PSFA.tif'));
    [PSFB, ~] = ImageJ_formatted_TIFF.ReadTifStack(strcat(path_psf, 'PSFB.tif'));
    [PSFA_bp, ~] = ImageJ_formatted_TIFF.ReadTifStack(strcat(path_psf, 'PSFA_BP.tif'));
    [PSFB_bp, ~] = ImageJ_formatted_TIFF.ReadTifStack(strcat(path_psf, 'PSFB_BP.tif'));
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

    switch(bgChoice)
                % affine registration method: only if regChoice == 2, 3, 4
                % 0: no registration, but perform affine trans;
                % 1: translation only; 
                % 2: rigid body; 
                % 3: 7 degrees of freedom (translation, rotation, scaling equally in 3 dimemsions)
                % 4: 9 degrees of freedom (translation, rotation, scaling); 
                % 5: 12 degrees of freedom; 
                % 6: rigid body first, then 12 degrees of freedom
                % 7: translation, rigid body, 9 degrees of freedom and then 12 degrees of freedom
        case 'translation only'
            affMethod = 1;
        case 'rigid body'
            affMethod = 2;
        case '7 DOF'
            affMethod = 3;
        case '9 DOF'
            affMethod = 4;
        case '12 DOF'
            affMethod = 7; % 3 DOF --> 6 DOF --> 9 DOF --> 12 DOF
    end
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
        filename_A = strcat(path_data, 'StackA_', num2str(imgNum), '.tif');
        filename_B = strcat(path_data, 'StackB_', num2str(imgNum), '.tif');
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
end
