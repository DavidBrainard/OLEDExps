function LaunchTomeMappingDemo()
%  Demo illustrating xyY-based tone-mapping pipeline

    rootDir = fileparts(which(mfilename));
    cd(rootDir);
    addpath(genpath(rootDir));
    addpath(genpath('MatlabEXR'));
    
    tonemapper = ToneMapper();
end







