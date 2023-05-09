function lib_unload(libName)
%{
unload dynamic-link library from MATLAB
Nov. 20, 2020, Min Guo

****************** Input ******************
libName: the library name (without .dll, .so)
%}

% check if dynamic-link library is already in MATLAB
if(libisloaded(libName))
    unloadlibrary(libName);
else
    disp(strcat('lib_unload: library not in MATLAB, no unloading for: ', libName));
end