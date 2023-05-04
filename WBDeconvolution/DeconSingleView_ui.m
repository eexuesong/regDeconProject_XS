% Xuesong Li 04/17/2023:

%{
SingleViewDcon_ui.m (GUI version):
single-view RL deconvolution for time-lapse images, compatible with unmatched back projectors.

Suitable microscope types:
* widefield fluorescence microscopy;
* confocal microscopy;
* instant structured illumination microscopy (iSIM);
* light-sheet fluorescence microscopy;

Before running the code, users should organize all input data as TIFF stacks (either 16-bit or 32-bit) as shown in Fig. C1.
All time-series (4D) data are grouped into a single folder:
Stack_0.tif, Stack_1.tiff, Stack_2.tif, ..., etc,
where the number indicates the relevant time point.
%}

%{
Richardson–Lucy deconvolution (Fourier domain):
(1) Spatial domain: Blurred_estimate = conv(Estimate(t), PSF)
    Fourier domain: fft(Blurred_estimate) = fft(Estimate(t)) .* OTF = fft(Estimate(t)) .* fft(PSF)
Note: Estimate(t) means current estimate.

(2) Ratio = Observed_image ./ Blurred_estimate
Note: Observed_image is input_image here.

(3) Spatial domain: Correction = conv(Ratio, PSF_flip)
    Fourier domain: fft(Correction) = fft(Ratio) .* OTF_BP = fft(Ratio) .* fft(PSF_flip)

(4) Estimate(t + 1) = Estimate(t) .* Correction
Note: Estimate(t + 1) means output estimate.
%}

clear all;

%% GUI part
% load raw data
[filename_data, path_data] = uigetfile('*.tif','Choose any one of raw data');
% load PSF data
[filename_psf, path_psf] = uigetfile('*.tif','Choose any one of PSF image', path_data);

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
timepoints = strsplit(answer{4},'-');
t1 = str2double(timepoints{1});
t2 = str2double(timepoints{2});
gpuFlag = 0;
if proMode == 1
    gpuFlag = 1;    % 0: CPU; 1: GPU
    gpuDevNum = 1;  % specify the GPU device if there are multiple GPUs
end

%% Create an output folder
path_output = strcat(path_data, 'results\');
mkdir(path_output);

%% Read images
% stackIn = single(ReadTifStack([path_data, filename_data]));
[stack_In, header_data] = ImageJ_formatted_TIFF.ReadTifStack(strcat(path_data, filename_data));
[Sy, Sx, Sz] = size(stack_In);

%% Forward and back projectors
disp('Preprocessing forward and back projectors ...');

% forward projector: PSF
[PSF_In, header_psf] = ImageJ_formatted_TIFF.ReadTifStack(strcat(path_psf, filename_psf));
PSF_In = single(PSF_In);
PSF1 = PSF_In / sum(PSF_In(:));

% back projector: PSF_bp
% parameters: light sheet microscopy as an example
switch(deconMethod)
    case 1
        bp_type = 'traditional';
    case 2
        bp_type = 'wiener-butterworth';
    otherwise
        error('Processing terminated, please set deconconvolution method as 1 or 2')
end

alpha = 0.05;   % Wiener constant w^2 in Wiener filter. If 1: use OTF value of PSF_bp at resolution limit.
beta = 1;       % Beta value in Butterworth filter. If 1: use OTF value of PSF_bp at resolution limit.
n = 10;         % 10: order of the Butterworth filter
resFlag = 1;    % 0: use PSF_fp FWHM / root(2) as resolution limit (for iSIM);
% 1: use PSF_fp FWHM as resoltuion limit;
% 2: use input values (iRes) as resoltuion limit;
iRes = [2.44, 2.44, 10];    % input resolution limit in 3 dimensions in terms of pixels;
verboseFlag = 0;
[PSF2, ~] = BackProjector(PSF1, bp_type, alpha, beta, n, resFlag, iRes, verboseFlag);
PSF2 = single(PSF2);
PSF2 = PSF2 / sum(PSF2(:));

if isempty(header_psf.resolution)
    ImageJ_formatted_TIFF.WriteTifStack(PSF1, strcat(path_output, 'PSF_fp.tif'));
    ImageJ_formatted_TIFF.WriteTifStack(PSF2, strcat(path_output, 'PSF_bp.tif'));
elseif isempty(header_psf.spacing)
    switch header_psf.unit
        case 'um'
            ImageJ_formatted_TIFF.WriteTifStack(PSF1, strcat(path_output, 'PSF_fp.tif'), header_psf.resolution);
            ImageJ_formatted_TIFF.WriteTifStack(PSF1, strcat(path_output, 'PSF_bp.tif'), header_psf.resolution);
        case 'nm'
            ImageJ_formatted_TIFF.WriteTifStack(PSF1, strcat(path_output, 'PSF_fp.tif'), header_psf.resolution / 1000);
            ImageJ_formatted_TIFF.WriteTifStack(PSF1, strcat(path_output, 'PSF_bp.tif'), header_psf.resolution / 1000);
        otherwise
            ImageJ_formatted_TIFF.WriteTifStack(PSF1, strcat(path_output, 'PSF_fp.tif'), header_psf.resolution);
            ImageJ_formatted_TIFF.WriteTifStack(PSF1, strcat(path_output, 'PSF_bp.tif'), header_psf.resolution);
    end
else
    switch header_psf.unit
        case 'um'
            ImageJ_formatted_TIFF.WriteTifStack(PSF1, strcat(path_output, 'PSF_fp.tif'), header_psf.resolution, header_psf.spacing);
            ImageJ_formatted_TIFF.WriteTifStack(PSF2, strcat(path_output, 'PSF_bp.tif'), header_psf.resolution, header_psf.spacing);
        case 'nm'
            ImageJ_formatted_TIFF.WriteTifStack(PSF1, strcat(path_output, 'PSF_fp.tif'), header_psf.resolution / 1000, header_psf.spacing / 1000);
            ImageJ_formatted_TIFF.WriteTifStack(PSF2, strcat(path_output, 'PSF_bp.tif'), header_psf.resolution / 1000, header_psf.spacing / 1000);
        otherwise
            ImageJ_formatted_TIFF.WriteTifStack(PSF1, strcat(path_output, 'PSF_fp.tif'), header_psf.resolution, header_psf.spacing);
            ImageJ_formatted_TIFF.WriteTifStack(PSF2, strcat(path_output, 'PSF_bp.tif'), header_psf.resolution, header_psf.spacing);
    end
end

%% PSF padding to be consistent with data
PSF_fp = align_size(PSF1, Sx, Sy, Sz);
PSF_bp = align_size(PSF2, Sx, Sy, Sz);
if(gpuFlag)
    g = gpuDevice(gpuDevNum);
    reset(g); wait(g);
    disp(append('GPU Memory before RL deconvolution: ', num2str(g.FreeMemory / 1024 / 1024 / 1024), ' GB'));

    OTF_fp = fftn(ifftshift(gpuArray(single(PSF_fp))));
    OTF_bp = fftn(ifftshift(gpuArray(single(PSF_bp))));
else
    OTF_fp = fftn(ifftshift(PSF_fp));
    OTF_bp = fftn(ifftshift(PSF_bp));
end

%% Deconvolution
% set initialization of the deconvolution
constInitial_flag = 0;   % 0: input image; 1: constant mean

disp('Start deconvolution...');
% smallValue = 0.001;
for imgNum = t1:t2
    disp(append('...Processing image #: ', num2str(imgNum)));
    filename_In = strcat(path_data, 'Stack_', num2str(imgNum), '.tif');
    [stack_In, header_data] = ImageJ_formatted_TIFF.ReadTifStack(filename_In);
    stack_In = single(stack_In);

    if(gpuFlag)
        stack = gpuArray(stack_In);
    else
        stack = stack_In;
    end
    %     stack = max(stack, smallValue);
    stack(stack <= 0) = eps;

    if constInitial_flag == 1
        Estimate = ones([Sy, Sx, Sz], 'single') * mean(stack(:));  % constant initialization
    else
        Estimate = stack;   % Measured image as initialization
    end

    for i = 1:itNum
        %         stackEstimate = stackEstimate .* ConvFFT3_S(stack ./ ConvFFT3_S(stackEstimate, OTF_fp), OTF_bp);
        Blur = real(ifftn(fftn(Estimate) .* OTF_fp));
        Ratio = stack ./ Blur;
        Correction = real(ifftn(fftn(Ratio) .* OTF_bp));
        Estimate = Estimate .* Correction;

        %         Estimate = max(Estimate, smallValue);
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
    filename_Out = strcat(path_output, 'Decon_test', num2str(imgNum), '.tif');
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

if(gpuFlag)
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