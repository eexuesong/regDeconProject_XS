% Generate back projector
clear all;

%% Forward projectors
path_psf = 'D:\Code\Matlab_Code\Code_from_Min\regDeconProject-master_XS\WBDeconvolution\DataForTest_lightsheet';
filename_psf = 'PSF.tif';
path_output = 'D:\Code\Matlab_Code\Code_from_Min\regDeconProject-master_XS\WBDeconvolution\DataForTest_lightsheet';

[psf_fp, header_psf] = ImageJ_formatted_TIFF.ReadTifStack(strcat(path_psf, '\', filename_psf));
[Sy, Sx, Sz] = size(psf_fp);

%% Create an output folder
path_output = strcat(path_output, '\');

if isequal(exist(path_output, 'dir'), 7)
    disp(append('output folder: ', path_output));
else
    mkdir(path_output);
    disp(append('output folder created: ', path_output));
end

%% back projector: PSF_bp
%%%%%%% ********** Parameters *************** %%%%%%
bp_type = 'wiener-butterworth';
alpha = 0.05;               % Wiener constant w^2 in Wiener filter. If 1: use OTF value of PSF_bp at resolution limit.
beta = 1;                   % Beta value in Butterworth filter. If 1: use OTF value of PSF_bp at resolution limit.
n = 10;                     % 10: order of the Butterworth filter
resFlag = 2;                % 0: use PSF_fp FWHM / root(2) as resolution limit (for iSIM);
                            % 1: use PSF_fp FWHM as resoltuion limit;
                            % 2: use input values (iRes) as resoltuion limit;
% iRes = zeros(1, 3);
iRes = [2.44, 2.44, 10];    % input resolution limit in 3 dimensions in terms of pixels;
verboseFlag = 1;

%% call function: BackProjector
[psf_bp, otf_bp] = BackProjector(psf_fp, bp_type, alpha, beta, n, resFlag, iRes, verboseFlag);

%% Save back projectors
otf_bp = fftshift(abs(otf_bp)) / abs(otf_bp(1));

if isempty(header_psf.resolution)
    ImageJ_formatted_TIFF.WriteTifStack(psf_bp, strcat(path_output, 'PSF_bp_', bp_type, '.tif'));
    ImageJ_formatted_TIFF.WriteTifStack(otf_bp, strcat(path_output, 'OTF_bp_', bp_type, '.tif'));
elseif isempty(header_psf.spacing)
    switch header_psf.unit
        case 'um'
            ImageJ_formatted_TIFF.WriteTifStack(psf_bp, strcat(path_output, 'PSF_bp_', bp_type, '.tif'), header_psf.resolution);
            ImageJ_formatted_TIFF.WriteTifStack(otf_bp, strcat(path_output, 'OTF_bp_', bp_type, '.tif'), header_psf.resolution);
        case 'nm'
            ImageJ_formatted_TIFF.WriteTifStack(psf_bp, strcat(path_output, 'PSF_bp_', bp_type, '.tif'), header_psf.resolution / 1000);
            ImageJ_formatted_TIFF.WriteTifStack(otf_bp, strcat(path_output, 'OTF_bp_', bp_type, '.tif'), header_psf.resolution / 1000);
        otherwise
            ImageJ_formatted_TIFF.WriteTifStack(psf_bp, strcat(path_output, 'PSF_bp_', bp_type, '.tif'), header_psf.resolution);
            ImageJ_formatted_TIFF.WriteTifStack(otf_bp, strcat(path_output, 'OTF_bp_', bp_type, '.tif'), header_psf.resolution);
    end
else
    switch header_psf.unit
        case 'um'
            ImageJ_formatted_TIFF.WriteTifStack(psf_bp, strcat(path_output, 'PSF_bp_', bp_type, '.tif'), header_psf.resolution, header_psf.spacing);
            ImageJ_formatted_TIFF.WriteTifStack(otf_bp, strcat(path_output, 'OTF_bp_', bp_type, '.tif'), header_psf.resolution, header_psf.spacing);
        case 'nm'
            ImageJ_formatted_TIFF.WriteTifStack(psf_bp, strcat(path_output, 'PSF_bp_', bp_type, '.tif'), header_psf.resolution / 1000, header_psf.spacing / 1000);
            ImageJ_formatted_TIFF.WriteTifStack(otf_bp, strcat(path_output, 'OTF_bp_', bp_type, '.tif'), header_psf.resolution / 1000, header_psf.spacing / 1000);
        otherwise
            ImageJ_formatted_TIFF.WriteTifStack(psf_bp, strcat(path_output, 'PSF_bp_', bp_type, '.tif'), header_psf.resolution, header_psf.spacing);
            ImageJ_formatted_TIFF.WriteTifStack(otf_bp, strcat(path_output, 'OTF_bp_', bp_type, '.tif'), header_psf.resolution, header_psf.spacing);
    end
end



