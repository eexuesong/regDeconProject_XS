function lib_load(libPath, libName, headFileName)
%{
load dynamic-link library to MATLAB cross-platforms: Windows, Linux, and Mac
Nov. 20, 2020, Min Guo

****************** Input ******************
libPath: path to the library
libName: the library name (without .dll, .so)
headFileName: the head file name (without .h), if not set, same as libName
%}

switch(nargin)
    case 2
        headFileName = libName;
    case 3
        % do nothing
    otherwise
        error('lib_load: unmatched arguments, please input 2 or 3 arguments!')
end

% check if dynamic-link library is already in MATLAB
if(libisloaded(libName))
    disp(strcat('lib_load: library already in MATLAB, no loading for: ', libName))
    return;
end

% load dynamic-link library
if ispc
    libFile = fullfile(libPath, strcat(libName, '.dll'));
elseif isunix
    libFile = fullfile(libPath, strcat(libName, '.so'));
elseif ismac
    libFile = fullfile(libPath, strcat(libName, '.so')); % ? ?
else
    error('lib_load: can not identify operation system');
end
libHFile = fullfile(libPath, strcat(headFileName, '.h'));
loadlibrary(libFile, libHFile);    