% Xuesong Li 05/08/2023:

%{
Registration.m: do registration for time-sequence images of two views,
the two view images are assumed to have same pixel size and oriented 
in the same direction.

Before running the code, users should organize all the input data as TIFF 
stacks (Fig. C1): StackA_N, StackB_N … etc, where the capital letters 
A and B indicate the two perpendicular views; and the number N indicates 
the relevant time point. Here we will register B view to A view.

While running the Registration.m script, two dialog boxes will pop up to 
guide users to 
(1) upload the raw dual-view data (see Fig. C1);
(2) set parameters (see Fig. C4) to 
a) choose registration mode:
    1: translation only;
    2: rigid body;
    3: 7 degrees of freedom (translation, rotation, scaling equally in 3 dimensions);
    4: 9 degrees of freedom (translation, rotation, scaling);
    5: 12 degrees of freedom (translation, rotation, scaling, shearing);
b) set the number of time points (in range form, e.g. 0-9) to be processed.

After registration, the registered data (i.e., B view) will be saved in the
same folder that contains the input data. They are saved as StackB_reg_N, 
where N indicates the number of time points.
%}

clear all;
warning('off', 'all');

%% GUI part
% load raw data
[filename_data, path_data] = uigetfile('*.tif', 'Choose any one of raw data');

% set parameters
inputDialog = uifigure('Name', 'Set Parameters', 'Position', [300 500 400 280]);
% set registration mode
uilabel(inputDialog, 'Text', 'Registration mode:', 'Position', [30 250 300 20]);
% Create a button group and radio buttons:
bg = uibuttongroup('Parent', inputDialog, 'Position', [50 130 300 120]);
rb1 = uiradiobutton(bg, 'Text', 'translation only', 'Position', [10 90 200 20]);
rb2 = uiradiobutton(bg, 'Text', 'rigid body', 'Position', [10 70 200 20]);
rb3 = uiradiobutton(bg, 'Text', '7 DOF', 'Position', [10 50 200 20]);
rb4 = uiradiobutton(bg, 'Text', '9 DOF', 'Position', [10 30 200 20]);
rb5 = uiradiobutton(bg, 'Text', '12 DOF', 'Position', [10 10 200 20]);
% for time points
uilabel(inputDialog, 'Text', 'Enter time points to be processed:', 'Position', [30 100 300 20]);
txa = uitextarea(inputDialog, 'Value', '0-1', 'Position', [50 70 100 30]);
% buttons to confirm/cancel processing
pb1 = uibutton(inputDialog, 'push', 'Position', [130 30 50 20], 'Text', 'Yes',...
    'ButtonPushedFcn', @(pb1, event) pbYes(inputDialog, bg, txa, path_data));
pb2 = uibutton(inputDialog, 'push', 'Position', [210 30 50 20], 'Text', 'Cancel',...
    'ButtonPushedFcn', @(pb1,event) pbCancel(inputDialog));

%% Cancel dialog box
function pbCancel(inputDialog)
    delete(inputDialog);
    disp('Processing cancelled!!!')
end

%% Yes dialog box
function pbYes(inputDialog, bg, txa, path_data)
    bgChoice = get(get(bg, 'SelectedObject'), 'Text');
    timeStr = get(txa, 'Value');
    timepoints = strsplit(timeStr{1}, '-');
    t1 = str2double(timepoints{1});
    t2 = str2double(timepoints{2});
    delete(inputDialog);    
    disp('Initialize processing...');

    %% Load DLL lib
    libPath = '..\cudaLib\';
    libName = 'libapi';

    %% Create an output folder
    path_output = strcat(path_data, 'results\');
    mkdir(path_output);

    %% Create arguments
    % configurations
    regChoice = 2; % *** registration choice: regChoice
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
            affMethod = 7;
    end
    flagTmx = 0;    % whether or not use input matrix, 1: yes; 0: no
                    % *** flagTmx: only if regChoice == 0, 2
    iTmx = eye(4);  % initial matrix for registration between two views. Tmx = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
    FTOL = 0.0001;  % reg threshold
    itLimit = 2000; % maximun iteration number for registration
    deviceNum = 0;  % GPU device: numbering from 0 by CUDA;
    %     verbose = 0;    % show details during processing: does not work for MATLAB
    %     records = zeros(1, 11);
    %     regRecords = libpointer('singlePtr', records);  % reg records and feedback
    gpuMemMode = 1; % 1: efficient GPU mode; 2: GPU memory-saved mode
    tic;
    
    %% Registration
    disp('Start registration...');
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

        % run registration function: reg3d
        [stackB_reg, oTmx] = reg3d_CUDA(stackA, stackB, libPath, libName, regChoice, affMethod, flagTmx, iTmx, FTOL, itLimit, deviceNum, gpuMemMode);
        
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
    cTime3 = toc;

    disp(append('... Total time cost: ', num2str(cTime3), ' s'));
    disp('Registration completed !!!');
end
