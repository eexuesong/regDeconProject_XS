% Xuesong Li 05/04/2023:

%{
DualViewDcon.m: dual-view joint RL deconvolution for time-lapse images, 
compatible with unmatched back projectors.

Suitable microscope types:
* dual-view light-sheet fluorescence microscopy (diSPIM);

For multiview data, users should organize the data as:
StackA_N, StackB_N, StackC_N, StackD_N, … etc, 
where the capital letters A/B/C indicate the relevant view and 
the number N indicates the relevant time point.

All views must be registered and have the same dimensions.

Users need to provide forward projectors (PSFs) of the microscope.
These PSFs should be oriented from the same perspective as their corresponding views.
They can be placed in the same folder as the input data or grouped into another folder,
named PSFA for the first view, PSFB for the second view, and so on (Fig. C2).

In the multiview deconvolution cases, the corresponding Wiener-Butterworth 
back projectors can be generated by Generate_BackProjector.m 
and fed to the deconvolution scripts correspondingly.
%}

clear all;

%% GUI part
% load raw data
[filename_data, path_data] = uigetfile('*.tif', 'Choose any one of raw data');
% load PSF data 
[filename_psf, path_psf] = uigetfile('*.tif', 'Choose any one of PSF image', path_data);

% set parameters
dlg_title = 'Set Parameters';
prompt = {'Enter deconvolution method: 1 for traditional decon; 2 for WB',...
    'Enter processing mode: 0 for CPU; 1 for GPU', 'Enter iteration number: ', ...
    'Enter time points to be processed'};
num_lines = 2;
defaultans = {'2', '1', '1', '0-2'};
answer = inputdlg(prompt, dlg_title, num_lines, defaultans);

% get parameters
deconMethod = str2double(answer{1});    % Deconvolution method: 1 for traditional deconvolution; 2 for Wiener-Butterworth deconvolution
proMode = str2double(answer{2});        % processing mode: 0 for CPU; 1 for GPU
itNum = str2double(answer{3});          % iteration number
timepoints = strsplit(answer{4}, '-');
t1 = str2double(timepoints{1});
t2 = str2double(timepoints{2});
gpuFlag = 0;
if proMode == 1
    gpuFlag = 1;    % 0: CPU; 1: GPU  
    gpuDevNum = 1;  % specify the GPU device if there are multiple GPUs
end

%% Create an output folder
path_output = strcat(path_data, 'results_diSPIM\');
mkdir(path_output);

%% Read images
% stackIn = single(ReadTifStack([path_data, filename_data]));
[stack_In, header_data] = ImageJ_formatted_TIFF.ReadTifStack(strcat(path_data, filename_data));
[Sy, Sx, Sz] = size(stack_In);

%% Forward and back projectors
disp('Preprocessing forward and back projectors ...');

% forward projectors: PSFA and PSFB
[PSF_In, ~] = ImageJ_formatted_TIFF.ReadTifStack(strcat(path_psf, 'PSFA.tif'));
PSF_In = single(PSF_In);
PSF1 = PSF_In / sum(PSF_In(:));
[PSF_In, ~] = ImageJ_formatted_TIFF.ReadTifStack(strcat(path_psf, 'PSFB.tif'));
PSF_In = single(PSF_In);
PSF2 = PSF_In / sum(PSF_In(:));

% back projector: PSF_bp
% parameters: light sheet microscopy as an example
switch(deconMethod)
    case 1
        PSF3 = flipPSF(PSF1);
        PSF4 = flipPSF(PSF2);
    case 2
        [PSF_In, ~] = ImageJ_formatted_TIFF.ReadTifStack(strcat(path_psf, 'PSFA_BP.tif'));
        PSF_In = single(PSF_In);
        PSF3 = PSF_In / sum(PSF_In(:));
        [PSF_In, ~] = ImageJ_formatted_TIFF.ReadTifStack(strcat(path_psf, 'PSFB_BP.tif'));
        PSF_In = single(PSF_In);
        PSF4 = PSF_In / sum(PSF_In(:));
    otherwise
        error('Processing terminated, please set deconconvolution method as 1 or 2')
end

%% PSF padding to be consistent with data
PSFA_fp = align_size(PSF1, Sx, Sy, Sz);
PSFB_fp = align_size(PSF2, Sx, Sy, Sz);
PSFA_bp = align_size(PSF3, Sx, Sy, Sz);
PSFB_bp = align_size(PSF4, Sx, Sy, Sz);
if gpuFlag
    g = gpuDevice(gpuDevNum); 
    reset(g); wait(g); 
    disp(append('GPU Memory before RL deconvolution: ', num2str(g.FreeMemory / 1024 / 1024 / 1024), ' GB'));

    OTFA_fp = fftn(ifftshift(gpuArray(single(PSFA_fp))));
    OTFB_fp = fftn(ifftshift(gpuArray(single(PSFB_fp))));
    OTFA_bp = fftn(ifftshift(gpuArray(single(PSFA_bp))));
    OTFB_bp = fftn(ifftshift(gpuArray(single(PSFB_bp))));
else
    OTFA_fp = fftn(ifftshift(PSFA_fp));
    OTFB_fp = fftn(ifftshift(PSFB_fp));
    OTFA_bp = fftn(ifftshift(PSFA_bp));
    OTFB_bp = fftn(ifftshift(PSFB_bp));
end

%% Deconvolution
disp('Start deconvolution...');
% smallValue = 0.01;
for imgNum = t1:t2
    disp(append('...Processing image #: ', num2str(imgNum)));
    filename_In_A = strcat(path_data, 'StackA_', num2str(imgNum), '.tif');
    filename_In_B = strcat(path_data, 'StackB_', num2str(imgNum), '.tif');
    [stack_In_A, header_data_A] = ImageJ_formatted_TIFF.ReadTifStack(filename_In_A);
    [stack_In_B, header_data_B] = ImageJ_formatted_TIFF.ReadTifStack(filename_In_B);
    stack_In_A = single(stack_In_A);
    stack_In_B = single(stack_In_B);

    if gpuFlag
        stackA = gpuArray(stack_In_A);
        stackB = gpuArray(stack_In_B);
    else
        stackA = stack_In_A;
        stackB = stack_In_B;
    end
    %     stackA = max(stackA, smallValue);
    %     stackB = max(stackB, smallValue);
    stackA(stackA <= 0) = eps;
    stackB(stackB <= 0) = eps;

    Estimate = (stackA + stackB) / 2;
    for i = 1:itNum
        %         Estimate = Estimate .* ConvFFT3_S(stackA ./ ConvFFT3_S(Estimate, OTFA_fp), OTFA_bp);
        %         Estimate = max(Estimate, smallValue);
        %         Estimate = Estimate .* ConvFFT3_S(stackB ./ ConvFFT3_S(Estimate, OTFB_fp), OTFB_bp);
        %         Estimate = max(Estimate, smallValue);
        
        % View A
        Blur = real(ifftn(fftn(Estimate) .* OTFA_fp));
        Ratio = stackA ./ Blur;
        Correction = real(ifftn(fftn(Ratio) .* OTFA_bp));
        Estimate = Estimate .* Correction;
        Estimate(Estimate <= 0) = eps;

        % View B
        Blur = real(ifftn(fftn(Estimate) .* OTFB_fp));
        Ratio = stackB ./ Blur;
        Correction = real(ifftn(fftn(Ratio) .* OTFB_bp));
        Estimate = Estimate .* Correction;
        Estimate(Estimate <= 0) = eps;
    end

    if gpuFlag
        stack_out = gather(Estimate);
    else
        stack_out = Estimate;
    end

    if header_data.BitsPerSample == 16
        if max(stack_out, [], 'all') <= 65535
            stack_out = uint16(stack_out);
        else
            stack_out = uint16(65535 * stack_out ./ max(stack_out, [], 'all'));
        end
    end
    
    % Write deconvolved images
    filename_Out = strcat(path_output, 'Decon_', num2str(imgNum), '.tif');
    if isempty(header_data.resolution)
        ImageJ_formatted_TIFF.WriteTifStack(stack_out, filename_Out);
    elseif isempty(header_data.spacing)
        switch header_data.unit
            case 'um'
                ImageJ_formatted_TIFF.WriteTifStack(stack_out, filename_Out, header_data.resolution);
            case 'nm'
                ImageJ_formatted_TIFF.WriteTifStack(stack_out, filename_Out, header_data.resolution / 1000);
            otherwise
                ImageJ_formatted_TIFF.WriteTifStack(stack_out, filename_Out, header_data.resolution);
        end
    else
        switch header_data.unit
            case 'um'
                ImageJ_formatted_TIFF.WriteTifStack(stack_out, filename_Out, header_data.resolution, header_data.spacing);
            case 'nm'
                ImageJ_formatted_TIFF.WriteTifStack(stack_out, filename_Out, header_data.resolution / 1000, header_data.spacing / 1000);
            otherwise
                ImageJ_formatted_TIFF.WriteTifStack(stack_out, filename_Out, header_data.resolution, header_data.spacing);
        end
    end
end

if gpuFlag
    reset(g);   % reset GPU
end

disp('Deconvolution completed !!!');

%% Function: zero padding PSF
function psf_padding = align_size(psf, nx_image, ny_image, nz_image, padValue)
if(nargin == 4)
    padValue = 0;
end

[ny_psf, nx_psf, nz_psf] = size(psf);
nx_max = max(nx_psf, nx_image);
ny_max = max(ny_psf, ny_image);
nz_max = max(nz_psf, nz_image);
imgTemp = ones(ny_max, nx_max, nz_max) * padValue;

nx_shift = round((nx_max - nx_psf) / 2) + 1;
ny_shift = round((ny_max - ny_psf) / 2) + 1;
nz_shift = round((nz_max - nz_psf) / 2) + 1;
imgTemp(ny_shift:ny_shift + ny_psf - 1, nx_shift:nx_shift + nx_psf - 1, nz_shift:nz_shift + nz_psf - 1) = psf;

nx_shift = round((nx_max - nx_image) / 2) + 1;
ny_shift = round((ny_max - ny_image) / 2) + 1;
nz_shift = round((nz_max - nz_image) / 2) + 1;
psf_padding = imgTemp(ny_shift:ny_shift + ny_image - 1, nx_shift:nx_shift + nx_image - 1, nz_shift:nz_shift + nz_image - 1);
end

%% Function: flip PSF relative to image center
function PSF_out = flipPSF(PSF_in)
% function outPSF = flipPSF(inPSF)
% outPSF(i, j, k) = inPSF(m - i + 1, n - j + 1, l - k + 1);
%{
Deprecated: Min's code
[ny, nx, nz] = size(PSF_in);
PSF_out = zeros(ny, nx, nz);
for i = 1:ny
    for j = 1:nx
        for k = 1:nz
            PSF_out(i, j, k) = PSF_in(ny - i + 1, nx - j + 1, nz - k + 1);
        end
    end
end
%}

PSF_out = flip(PSF_in, 1);
PSF_out = flip(PSF_out, 2);
PSF_out = flip(PSF_out, 3);
end