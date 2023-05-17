function [PSF_bp, OTF_bp] = BackProjector(PSF_fp, bp_type, alpha, beta, n, resFlag, iRes, verboseFlag)
% Generate Backprojector: PSF and OTF
% April 12, 2019 (Min Guo)
% Modified by Xuesong Li 04/17/2023

% % % Output
% PSF_bp: Back projector in spatial domain
% OTF_bp: Back projector in Fourier domain

% % % Input
% PSF_fp: Forward projector
% bp_type: 'traditional', 'gaussian', 'butterworth', 'wiener', 'wiener-butterworth'
% alpha: 0.0001 ~ 0.001
%      1: use OTF value of PSF_bp at resolution limit;
%      else: alpha = input alpha;
% beta: 0.001 ~ 0.01
%      1: use OTF value of PSF_bp at resolution limit;
%      else: beta = input beta;
% n: 4 ~ 15
%      order of the Butterworth filter
% resFlag:
%      0: use PSF_fp FWHM/root(2) as resolution limit (for iSIM);
%      1: use PSF_fp FWHM as resoltuion limit;
%      2: use input values (iRes) as resoltuion limit;
% iRes: 1 x 3 array
%      input resolution limit in 3 dimensions in terms of pixels;
% verboseFlag:
%      0: hide log and intermediate results
%      1: show log and intermediate results

% Notes by Xuesong Li 04/17/2023
%{
Wiener Filter:
The transfer function of Wiener filter is defined as:
H(u, v) = O*(u, v) / (abs(O(u, v))^2 + w^2)

where:
O(u, v) is the optical transfer function
O*(u, v) is the conjugate of O(u, v)
w is the Wiener constant
%}

%{
Butterworth Lowpass Filter (BLPF):
In the field of Image Processing, Butterworth Lowpass Filter (BLPF) is used for image smoothing in the frequency domain.
It removes high-frequency noise from a digital image and preserves low-frequency components.
The transfer function of BLPF of order n is defined as:
H(u, v) = 1 / sqrt(1 + (D(u, v) / D0)^(2n))

where:
D0 is a positive constant. BLPF passes all the frequencies less than D0 value without attenuation and cuts off all the frequencies greater than it.
This D0 is the transition point between H(u, v) = 1 and H(u, v) = 0, so this is termed as cutoff frequency.
But instead of making a sharp cut-off (like, Ideal Lowpass Filter (ILPF)), it introduces a smooth transition from 1 to 0 to reduce ringing artifacts.

D(u, v) is the Euclidean Distance from any point (u, v) to the origin of the frequency plane, i.e.,
D(u, v) = sqrt(u^2 + v^2)

n is the order of the Butterworth filter.
%}

%{
More general format of frequency response for a Butterworth Filter:
H(u, v) = 1 / sqrt(1 + ε^2 * (D(u, v) / D0)^(2n))

where:
n is the filter order.
Epsilon ε is the maximum pass band gain Amax. 
If Amax is defined at a frequency equal to the cut-off -3dB corner point (ƒc), ε will then be equal to one and therefore ε^2 will also be one.
However, if you now wish to define Amax at a different gain value, for example 1dB, or 1.1220 (1dB = 20 * log10(Amax)) then the new value of epsilon, ε is found by:
    H1 = H0 / sqrt(1 + ε^2)     -->     ε^2 = (H0 / H1)^2 - 1 = 1.122 ^ 2 - 1 = 0.2589
where:
H0 = the Maximum Pass band Gain, Amax
H1 = the Minimum Pass band Gain.

If our code for 'butterworth', H0 / H1 = 1 / beta. So, ε^2 = 1 / beta^2 - 1;
Then H(u, v) = 1 / sqrt(1 + (1 / beta^2 - 1) * (D(u, v) / D0)^(2n))
When D(u, v) = D0, the cutoff frequency, H(u, v) = beta.

If our code for 'wiener-butterworth', H0 / H1 = beta_wiener_x / beta. So, ε^2 = (beta_wiener_x / beta)^2 - 1;
Then H(u, v) = 1 / sqrt(1 + ((beta_wiener_x / beta)^2 - 1) * (D(u, v) / D0)^(2n))
When D(u, v) = D0, the cutoff frequency, H(u, v) = beta / beta_wiener_x.
%}



if(nargin == 1)
    bp_type = 'traditional';
    alpha = 0.001;
    beta = 1;
    n = 10;
    resFlag = 1;
    iRes = [0,0,0];
    verboseFlag = 0;
end
if(nargin == 2)
    alpha = 0.001;
    beta = 1;
    n = 10;
    resFlag = 1;
    iRes = [0,0,0];
    verboseFlag = 0;
end
if(nargin == 3)
    beta = 1;
    n = 10;
    resFlag = 1;
    iRes = [0,0,0];
    verboseFlag = 0;
end
if(nargin == 4)
    n = 10;
    resFlag = 1;
    iRes = [0,0,0];
    verboseFlag = 0;
end
if(nargin == 5)
    resFlag = 1;
    iRes = [0,0,0];
    verboseFlag = 0;
end
if(nargin == 7)
    verboseFlag = 0;
end

%% Input PSF size and center
PSF_fp = single(PSF_fp);
[ny, nx, nz] = size(PSF_fp);

nxo = (nx + 1) / 2;
nyo = (ny + 1) / 2;
nzo = (nz + 1) / 2;

Sox = round(nxo);
Soy = round(nyo);
Soz = round(nzo);

if verboseFlag
    disp(append('Back projector type: ', bp_type));
end

%% Calculate PSF and OTF size
[FWHM_x, FWHM_y, FWHM_z] = fwhm_PSF(PSF_fp);
if(verboseFlag)
    disp(append('Forward projector FWHMs [nx, ny, nz]: ', num2str(FWHM_x), ' x ', num2str(FWHM_y), ' x ', num2str(FWHM_z)));
end

%% Set resolution cutoff
switch(resFlag)
    case 0  % Set resolution as 1 / root(2) of PSF_fp FWHM: iSIM case
        resx = FWHM_x / sqrt(2); resy = FWHM_y / sqrt(2); resz = FWHM_z / sqrt(2);
    case 1  % Set resolution as PSF_fp FWHM
        resx = FWHM_x; resy = FWHM_y; resz = FWHM_z;
    case 2  % Set resolution based input values
        resx = iRes(1); resy = iRes(2); resz = iRes(3);
    otherwise
        error('Processing terminated, please set resFlag as 0, 1, or 2')
end
% pixel size in Fourier domain
px = 1 / nx;
py = 1 / ny;
pz = 1 / nz;
% frequency cutoff in terms of pixels (radius)
OTF_cutoff_x = 1 / resx / px;
OTF_cutoff_y = 1 / resy / py;
OTF_cutoff_z = 1 / resz / pz;
if(verboseFlag)
    disp(['Resolution cutoff in spatial domain [nx, ny, nz]: ' num2str(resx) ' x ' num2str(resy) ' x ' num2str(resz)]);
    disp(['Resolution cutoff in Fourier domain [kx, ky, kz]: ' num2str(OTF_cutoff_x) ' x ' num2str(OTF_cutoff_y) ' x ' num2str(OTF_cutoff_z)]);
end

%% Normalize OTF: traditional back projector
PSF_flip = flipPSF(PSF_fp);
OTF_flip = fftn(ifftshift(PSF_flip));
OTF_flip_abs = fftshift(abs(OTF_flip));
OTF_flip_max = max(OTF_flip_abs(:));    % find maximum value
M = OTF_flip_max(1);
OTF_flip_abs_norm = OTF_flip_abs / M;

%% Check cutoff gains of traditional back projector
OTF_plane_xy_max = squeeze(max(OTF_flip_abs_norm, [], 3));
OTF_line_x_max = max(OTF_plane_xy_max, [], 1);
to1 = max(round(nxo - OTF_cutoff_x), 1); to2 = min(round(nxo + OTF_cutoff_x), nx);
beta_fp_x = (OTF_line_x_max(to1) + OTF_line_x_max(to2)) / 2; % OTF frequency intensity at cutoff: x

OTF_plane_xy_max = squeeze(max(OTF_flip_abs_norm, [], 3));
OTF_line_x_max = max(OTF_plane_xy_max, [], 2);
to1 = max(round(nyo - OTF_cutoff_y), 1); to2 = min(round(nyo + OTF_cutoff_y), ny);
beta_fp_y = (OTF_line_x_max(to1) + OTF_line_x_max(to2)) / 2; % OTF frequency intensity at cutoff: y

OTF_plane_xz_max = squeeze(max(OTF_flip_abs_norm, [], 1));
OTF_line_z_max = max(OTF_plane_xz_max, [], 1);
to1 = max(round(nzo - OTF_cutoff_z), 1); to2 = min(round(nzo + OTF_cutoff_z), nz);
beta_fp_z = (OTF_line_z_max(to1) + OTF_line_z_max(to2)) / 2; % OTF frequency intensity at cutoff: z

beta_fp = (beta_fp_x + beta_fp_y + beta_fp_z) / 3;
if(verboseFlag)
    disp(['Cutoff gain of forward projector: ' num2str(beta_fp_x) ' x ' num2str(beta_fp_y)...
        ' x ' num2str(beta_fp_z) ', Average = ' num2str(beta_fp)]);
end

%% Parameter for Wiener filter
if (alpha == 1)
    alpha = beta_fp;
    if(verboseFlag)
        disp(['Wiener parameter adjusted as traditional BP cutoff gain: alpha = ' num2str(alpha)]);
    end
else
    if(verboseFlag)
        disp(['Wiener parameter set as input: alpha = ' num2str(alpha)]);
    end
end

%% Parameter for Butterworth filter
if(beta == 1)
    beta = beta_fp;
    if(verboseFlag)
        disp(['Cutoff gain adjusted as traditional BP cutoff gain: beta = ' num2str(beta)]);
    end
else
    if(verboseFlag)
        disp(['Cutoff gain set as input: beta = ' num2str(beta)]);
    end
end

%% Order of Butterworth filter
% pn = 2 * n;
if(verboseFlag)
    disp(['Butterworth order (slope parameter) set as: n = ' num2str(n)]);
end

%% Core code
switch bp_type
    case 'traditional'
        PSF_bp= flipPSF(PSF_fp);
        OTF_bp = fftn(ifftshift(PSF_bp));

    case 'gaussian'
        resx = FWHM_x;
        resy = FWHM_y;
        resz = FWHM_z;
        PSF_bp= gen_gaussianPSF_3D(nx, ny, nz, resx, resy, resz);
        OTF_bp = fftn(ifftshift(PSF_bp));

    case 'butterworth'
        %{
        OTF_butterworth = 1 / sqrt(1 + ee * (kx / kcx)^pn)
        beta = 1 / sqrt(1 + ee) --> ee = 1 / beta^2 - 1;
        %}
        % create Butteworth Filter
        kcx = OTF_cutoff_x; % width of Butterworth Filter
        kcy = OTF_cutoff_y; % width of Butterworth Filter
        kcz = OTF_cutoff_z; % width of Butterworth Filter
        ee = 1 / beta^2 - 1;
        %{
        Deprecated: Min's code
        mask = zeros(ny, nx, nz);
        for i = 1: ny
            for j = 1 : nx
                for k = 1 : nz
                    w = ((i - nyo) / kcy)^2 + ((j - nxo) / kcx)^2 + ((k - nzo) / kcz)^2;
                    mask(i, j, k) = 1 / sqrt(1 + ee * w^n); % w^n = (kx / kcx)^pn
                end
            end
        end
        mask = 1 / sqrt(1 + ee * w^n);
        %}

        x = 1 : nx;
        y = 1 : ny;
        z = 1 : nz;
        [X, Y, Z] = meshgrid(x, y, z);
        W = ((X - nxo) / kcx).^2 + ((Y - nyo) / kcy).^2 + ((Z - nzo) / kcz).^2;
        Mask = 1 / sqrt(1 + ee * W.^n);

        OTF_bp = ifftshift(Mask);
        PSF_bp = fftshift(real(ifftn(OTF_bp)));

    case 'wiener'
        OTF_flip_norm = OTF_flip / M;   % Normalized OTF_flip
        OTF_bp = OTF_flip_norm ./ (abs(OTF_flip_norm).^2 + alpha);   % Wiener filter
        PSF_bp = fftshift(real(ifftn(OTF_bp)));

    case 'wiener-butterworth'
        % create Wiener filter
        OTF_flip_norm = OTF_flip / M;
        OTF_Wiener = OTF_flip_norm ./ (abs(OTF_flip_norm).^2 + alpha);
        % get cutoff gain of Wiener filter
        OTF_Wiener_abs = fftshift(abs(OTF_Wiener));
        OTF_Wiener_abs_xy = abs(squeeze(OTF_Wiener_abs(:, :, Soz))); % central XY slice
        OTF_Wiener_abs_x_max = max(OTF_Wiener_abs_xy, [], 1);
        to1 = max(round(nxo - OTF_cutoff_x), 1); to2 = min(round(nxo + OTF_cutoff_x), nx);
        beta_wiener_x = (OTF_Wiener_abs_x_max(to1) + OTF_Wiener_abs_x_max(to2)) / 2;    % OTF frequency intensity at cutoff: x
        if(verboseFlag)
            disp(['Wiener cutoff gain: beta_wiener_x = ' num2str(beta_wiener_x)]);
        end

        %{
        OTF_wiener-butterworth = Wiener .* 1 / sqrt(1 + ee * (kx / kcx)^pn)
        beta = beta_wiener_x * 1 / sqrt(1 + ee) --> ee = (beta_wiener_x / beta)^2 - 1;
        %}
        % create Butteworth Filter
        kcx = OTF_cutoff_x; % width of Butterworth Filter
        kcy = OTF_cutoff_y; % width of Butterworth Filter
        kcz = OTF_cutoff_z; % width of Butterworth Filter
        ee = (beta_wiener_x / beta)^2 - 1;
        %{
        Deprecated: Min's code
        mask = zeros(ny, nx, nz);
        for i = 1 : ny
            for j = 1 : nx
                for k = 1 : nz
                    w = ((i - nyo) / kcy)^2 + ((j - nxo) / kcx)^2 + ((k - nzo) / kcz)^2;
                    mask(i, j, k) = 1 / sqrt(1 + ee * w^n); % w^n = (kx / kcx)^pn
                end
            end
        end
        mask = ifftshift(mask); % Butterworth Filter
        %}

        x = 1 : nx;
        y = 1 : ny;
        z = 1 : nz;
        [X, Y, Z] = meshgrid(x, y, z);
        W = ((X - nxo) / kcx).^2 + ((Y - nyo) / kcy).^2 + ((Z - nzo) / kcz).^2;
        Mask = 1 / sqrt(1 + ee * W.^n);
        Mask = ifftshift(Mask); % Butterworth Filter

        % create Wiener-Butteworth Filter
        OTF_bp = Mask .* OTF_Wiener;    % final OTF_bp cutfoff gain: beta
        PSF_bp = fftshift(real(ifftn(OTF_bp)));
    otherwise
        error('bp_type does not match any back-projector type')
end

if(verboseFlag)
    line1 = squeeze(PSF_flip(:, Sox, Soz));
    line2 = squeeze(PSF_bp(:, Sox, Soz));
    figure,subplot(1, 2, 1);
    plot(1 : ny, line1 / max(line1(:)), 1 : ny, line2 / max(line2(:)), 'LineWidth', 2);
    xlabel('Pixel Position');
    ylabel('Normalized Value');
    legend('Trad. bp', append(bp_type, ' bp'));
    title('Back Projector Profiles (Spatial Domain)');

    line1 = squeeze(abs(OTF_flip(1:Soy, 1, 1)));
    line2 = squeeze(abs(OTF_bp(1:Soy, 1, 1)));
    subplot(1, 2, 2);
    plot(0 : Soy - 1, line1 / max(line1(1)), 0 : Soy - 1, line2 / max(line2(1)), 'LineWidth', 2);
    hold on, plot(ones(1, 13) * OTF_cutoff_x, 0 : 0.1 : 1.2, '--r.');
    xlabel('Pixel Position');
    ylabel('Normalized Value');
    legend('Trad. bp', append(bp_type, ' bp'), 'resolution limit');
    title('Back Projectors Profiles (Frequency Domain)');
    disp('Back projector generated!!!');
end
end



%% Function: calculate FWHM of PSF
function [FWHM_x, FWHM_y, FWHM_z] = fwhm_PSF(PSF, pixelSize, cFlag, fitFlag)
% Feed back the full width at half maximun of the input PSF
% fwhm.m and mygaussfit.m are needed
% cFlag
%       0: use maximum's position as PSF center position
%       1: use matrix's center position as PSF center position
% fitFlag
%       0: no fitting before calculate FWHM
%       1: spine fitting before calculate FWHM
%       2: gaussian fitting before calculate FWHM
%
if(nargin == 1)
    pixelSize = 1;
    cFlag = 0;
    fitFlag = 0;
end

if(nargin == 2)
    cFlag = 0;
    fitFlag = 0;
end

if(nargin == 3)
    fitFlag = 0;
end

% PSF = PSF - mean(PSF(:));
[ny, nx, nz] = size(PSF);
if (ny == 1) || (nx == 1)
    % 1D input
    x = 1:max(ny, nx);
    x = x';
    y = PSF(:);
    if fitFlag == 1
        % fitFlag = 1: spine fitting before calculate FWHM
        xq = 1:0.1:nx;
        yq = interp1(x, y, xq, 'spline');
        FWHM_x = fwhm(xq, yq);
    elseif fitFlag == 2
        % fitFlag = 2: gaussian fitting before calculate FWHM
        [sig, ~, ~] = mygaussfit(x, y);
        %         FWHM_x = sig * 2.3548;
        FWHM_x = sig * 2 * sqrt(2 * log(2));
    else
        % fitFlag = 0: no fitting before calculate FWHM
        FWHM_x = fwhm(x, y);
    end

    FWHM_y = 0;
    FWHM_z = 0;
elseif nz == 1
    % 2D input
    if(cFlag)
        % cFlag = 1: use matrix's center position as PSF center position
        indy = floor((ny + 1) / 2);
        indx = floor((nx + 1) / 2);
    else
        % cFlag = 0: use maximum's position as PSF center position
        [~, ind] = max(PSF(:)); % find maximum value and position
        [indy, indx] = ind2sub([ny, nx], ind(1));
    end

    % Y dimension
    x = 1:ny;
    x = x';
    y = PSF(:, indx);
    y = y(:);
    if fitFlag == 1
        xq = 1:0.1:ny;
        yq = interp1(x, y, xq, 'spline');
        FWHM_y = fwhm(xq, yq);
    elseif fitFlag == 2
        [sig, ~, ~] = mygaussfit(x, y);
        %         FWHM_y = sig * 2.3548;
        FWHM_y = sig * 2 * sqrt(2 * log(2));
    else
        FWHM_y = fwhm(x, y);
    end

    % X dimension
    x = 1:nx;
    x = x';
    y = PSF(indy, :);
    y = y(:);
    if fitFlag == 1
        xq = 1:0.1:nx;
        yq = interp1(x, y, xq, 'spline');
        FWHM_x = fwhm(xq, yq);
    elseif fitFlag == 2
        [sig, ~, ~] = mygaussfit(x, y);
        %         FWHM_x = sig * 2.3548;
        FWHM_x = sig * 2 * sqrt(2 * log(2));
    else
        FWHM_x = fwhm(x, y);
    end

    FWHM_z = 0;
else
    % 3D input
    if(cFlag)
        indy = floor((ny + 1) / 2);
        indx = floor((nx + 1) / 2);
        indz = floor((nz + 1) / 2);
    else
        [~, ind] = max(PSF(:)); % find maximum value and position
        [indy, indx, indz] = ind2sub([ny, nx, nz], ind(1));
    end

    % Y dimension
    x = 1:ny;
    x = x';
    y = PSF(:, indx, indz);
    y = y(:);
    if fitFlag == 1
        xq = 1:0.1:ny;
        yq = interp1(x, y, xq, 'spline');
        FWHM_y = fwhm(xq, yq);
    elseif fitFlag == 2
        [sig, ~, ~] = mygaussfit(x, y);
        %         FWHM_y = sig * 2.3548;
        FWHM_y = sig * 2 * sqrt(2 * log(2));
    else
        FWHM_y = fwhm(x, y);
    end

    % X dimension
    x = 1:nx;
    x = x';
    y = PSF(indy, :, indz);
    y = y(:);
    if fitFlag == 1
        xq = 1:0.1:nx;
        yq = interp1(x, y, xq, 'spline');
        FWHM_x = fwhm(xq, yq);
    elseif fitFlag == 2
        [sig, ~, ~] = mygaussfit(x, y);
        %         FWHM_x = sig * 2.3548;
        FWHM_x = sig * 2 * sqrt(2 * log(2));
    else
        FWHM_x = fwhm(x, y);
    end

    % Z dimension
    x = 1:nz;
    x = x';
    y = PSF(indy, indx, :);
    y = y(:);
    if fitFlag == 1
        xq = 1:0.1:nz;
        yq = interp1(x, y, xq, 'spline');
        FWHM_z = fwhm(xq, yq);
    elseif fitFlag == 2
        [sig, ~, ~] = mygaussfit(x, y);
        %         FWHM_z = sig * 2.3548;
        FWHM_z = sig * 2 * sqrt(2 * log(2));
    else
        FWHM_z = fwhm(x, y);
    end
end

FWHM_x = FWHM_x * pixelSize;
FWHM_y = FWHM_y * pixelSize;
FWHM_z = FWHM_z * pixelSize;
end

%% Function: FWHM of the waveform y(x) and its polarity
function width = fwhm(x, y)
%{
Full-Width at Half-Maximum (FWHM) of the waveform y(x) and its polarity.
The FWHM result in 'width' will be in units of 'x'

Rev 1.2, April 2006 (Patrick Egan)
%}

y = y / max(y);
N = length(y);
lev50 = 0.5;

% find index of center (max or min) of pulse
if y(1) < lev50
    [~, centerindex] = max(y);
    Pol = +1;
    %     disp('Pulse Polarity = Positive')
else
    [~, centerindex] = min(y);
    Pol = -1;
    %     disp('Pulse Polarity = Negative')
end

% first crossing is between y(i-1) & y(i)
i = 2;
while sign(y(i) - lev50) == sign(y(i - 1) - lev50)
    i = i + 1;
end
interp = (lev50 - y(i - 1)) / (y(i) - y(i - 1));
tlead = x(i - 1) + interp * (x(i) - x(i - 1));

% start search for next crossing at center
i = centerindex + 1;
while ((sign(y(i) - lev50) == sign(y(i - 1) - lev50)) && (i <= N - 1))
    i = i + 1;
end

if i ~= N
    Ptype = 1;
    %     disp('Pulse is Impulse or Rectangular with 2 edges')
    interp = (lev50 - y(i - 1)) / (y(i) - y(i - 1));
    ttrail = x(i - 1) + interp * (x(i) - x(i - 1));
    width = ttrail - tlead;
else
    Ptype = 2;
    %     disp('Step-Like Pulse, no second edge')
    ttrail = NaN;
    width = NaN;
end
end

%% Function: 1-D Guaussian fitting
function [sigma, mu, A] = mygaussfit(x, y, h)
% [sigma,mu,A] = mygaussfit(x,y)
% [sigma,mu,A] = mygaussfit(x,y,h)

%{
This function is doing fit to the function
y = A * exp(-(x - mu)^2 / (2 * sigma^2))

The fitting is been done by a polyfit the lan of the data.

h is the threshold which is the fraction from the maximum y height that the data is been taken from.
h should be a number between 0-1.
If h have not been taken it is set to be 0.2 as default.
%}

% threshold
if nargin == 2
    h = 0.2;
end

% cutting
ymax = max(y);
xnew = [];
ynew = [];
for n = 1:length(x)
    if y(n) > ymax * h
        xnew = [xnew, x(n)];
        ynew = [ynew, y(n)];
    end
end

% fitting
ylog = log(ynew);
xlog = xnew;
p = polyfit(xlog, ylog, 2);
A2 = p(1);
A1 = p(2);
A0 = p(3);

%{
Note: ln(y) = -1 / (2 * sigma ^2) * x^2 + mu / sigma^2 * x + ln(A) - mu^2 / (2 * sigma^2)
% A2 = -1 / (2 * sigma ^2)
% A1 = mu / sigma^2
% A0 = ln(A) - mu^2 / (2 * sigma^2)
%}

sigma = sqrt(-1/(2 * A2));
mu = A1 * sigma^2;
A = exp(A0 + mu^2 / (2 * sigma^2));
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

%% Function: generate 3D Gaussian PSF with FWHM input
function PSF = gen_gaussianPSF_3D(nx, ny, nz, FWHM_x, FWHM_y, FWHM_z)
if(nargin == 4)
    FWHM_y = FWHM_x;
    FWHM_z = FWHM_x;
end
% sig_x = FWHM_x / 2.3548;
% sig_y = FWHM_y / 2.3548;
% sig_z = FWHM_z / 2.3548;
sig_x = FWHM_x / (2 * sqrt(2 * log(2)));
sig_y = FWHM_y / (2 * sqrt(2 * log(2)));
sig_z = FWHM_z / (2 * sqrt(2 * log(2)));
PSF = gen_gaussian3D(nx, ny, nz, sig_x, sig_y, sig_z);
end

%% Function: generate 3D Gaussian distribution with sig_x, sig_y, sig_z
function I = gen_gaussian3D(nx, ny, nz, sig_x, sig_y, sig_z)
%{
Multi-dimensional Gaussian function:
1 / (sigma_x * sqrt(2 * pi)) * exp(-(x - x0)^2 / (2 * sigma_x^2)) * ...
When it is a 3D guassian function:
f(x, y, z) = N * exp(-(x - x0)^2 / (2 * sigma_x^2) - (y - y0)^2 / (2 * sigma_y^2) - (z - z0)^2 / (2 * sigma_z^2))
where coefficient N = 1 / (sigma_x * sigma_y * sigma_z * (2 * pi)^(3 / 2))
%}
sqrSigx = sig_x^2 * 2;
sqrSigy = sig_y^2 * 2;
sqrSigz = sig_z^2 * 2;
nxo = (nx + 1) / 2;
nyo = (ny + 1) / 2;
nzo = (nz + 1) / 2;
coef = 1 / ((2 * pi)^(3 / 2) * sig_x * sig_y * sig_z);

%{
Deprecated: Min's code
I = zeros(ny, nx, nz, 'double');
for i = 1 : ny
    for j = 1 : nx
        for k = 1 : nz
            d = (i - nyo)^2 / sqrSigx + (j - nxo)^2 / sqrSigy + (k - nzo)^2 / sqrSigz;
            I(i, j, k) = exp(-d);
        end
    end
end
%}

x = 1 : nx;
y = 1 : ny;
z = 1 : nz;
[X, Y, Z] = meshgrid(x, y, z);
D = (X - nxo).^2 / sqrSigx + (Y - nyo).^2 / sqrSigy + (Z - nzo).^2 / sqrSigz;
I = exp(-D);
I = coef * I;
end


