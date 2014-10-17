function generateCloudProbe
    
    % Load the mean gamma curve function (gammaFunction)
    load('/Users/Shared/Matlab/Experiments/OLEDExps/PreliminaryData/GammaFunction.mat');
    
    stimParams.phaseVelocityAtZeroSF           = 0.5;      % phase velocity at zero SF
    stimParams.phaseVelocityFactor             = 0.5;      % velocity of phase at SF_i will be: SF_i ^ phaseVelocityFactor
    stimParams.meanChangeInPhaseAnglePerFrame  = pi/12;    % mean change in phase angle per frame
    stimParams.stdChangeInPhaseAnglePerFrame   = pi/18;
    stimParams.framesNum                       = 32;

    stimParams.deltaOri                 = 45;
    stimParams.exponentOfOneOverFArray  = [1 1.2 1.4];
    stimParams.oriBiasArray             = [2 1];
    stimParams.blockSizeArray           = [48 96 192];
    
    totalFrames = 0;
    for exponentOfOneOverFIndex = 1:numel(stimParams.exponentOfOneOverFArray)
        for oriBiasIndex = 1:numel(stimParams.oriBiasArray)
            [stimuli{exponentOfOneOverFIndex, oriBiasIndex}.imageSequence, ...
             stimuli{exponentOfOneOverFIndex, oriBiasIndex}.orientationSequence] = ...
                generateMovingCloudSequence(stimParams, ...
                stimParams.exponentOfOneOverFArray(exponentOfOneOverFIndex), ...
                stimParams.oriBiasArray(oriBiasIndex), ...
                gammaFunction);
            totalFrames = totalFrames + stimParams.framesNum;
        end
    end
    
    totalFrames
    
    % save stimuli
    save('PixelOLEDprobes.mat', 'stimParams', 'stimuli');
end

    
function [imageSequence, orientationSequence] = generateMovingCloudSequence(stimParams, exponentOfOneOverF, oribias, gammaFunction)
%    Create a movie of noise images having 1/f amplitude spectum properties
%
%    size       - size of image to produce
%    factor     - controls spectrum = 1/(f^factor)
%    meandev    - mean change in phaseangle per frame
%    stddev     - stddev of change in phaseangle per frame
%    lowvel     - phase velocity at 0 frequency 
%    velfactor  - phase velocity = freq^velfactor
%    nframes    - no of frames in movie
%
%    factor = 0             - raw Gaussian noise image
%           = 1             - gives the 1/f `standard' drop-off for `natural' images
%           = 1.5           - seems to give the most intersting `cloud patterns'
%           = 2 or greater  - produces `blobby' images

    % unload stimParams
    factor          = exponentOfOneOverF;
    lowvel          = stimParams.phaseVelocityAtZeroSF;
    velfactor       = stimParams.phaseVelocityFactor;
    meandev         = stimParams.meanChangeInPhaseAnglePerFrame;
    stddev          = stimParams.stdChangeInPhaseAnglePerFrame;
    nframes         = stimParams.framesNum;
    deltaOri        = stimParams.deltaOri;
    blockSizeArray  = stimParams.blockSizeArray;
   
    desiredMatrixSize = [1080 1920];
    filterSize  = 2048;
    
    ratio       = 2^oribias;
    orientation = 0;
    ellipse     = makeEllipseFilter(orientation, ratio, filterSize);
    filter      = 1./(ellipse.^factor);       % Construct the filter.
    filter      = fftshift(filter);
    
    phasemod    = fftshift(ellipse.^velfactor + lowvel);
    phase       = random('Uniform',0,2*pi, filterSize,filterSize);  % Random uniform distribution 0 - 2pi
   
    orientationSequence = [];
    if (oribias > 0)
        
        OrientationsTested = [0:deltaOri:180-deltaOri];
        NumberOfPresentationsForEachOrientation = 1;
        orientationSequence = randperm(length(OrientationsTested));
        for k = 1:NumberOfPresentationsForEachOrientation
            orientationSequence = [orientationSequence randperm(length(OrientationsTested))];
        end
        orientationSequence = (orientationSequence-1)*deltaOri;
        consecutiveFramesWithSimilarOrientation = nframes/(length(OrientationsTested)*NumberOfPresentationsForEachOrientation);
        orientationIndex = 0;
    end
    
    % pre-allocate memory
    imageSequence = zeros(2,nframes,desiredMatrixSize(1), desiredMatrixSize(2), 'uint8');
    
    for frameIndex = 1:nframes  
      if (oribias > 0)    
        if (mod((frameIndex-1),consecutiveFramesWithSimilarOrientation) == 0)      
            orientationIndex   = orientationIndex + 1;
            currentOrientation = orientationSequence(orientationIndex);
            % compute new filter
            ellipse = makeEllipseFilter(currentOrientation, oribias, filterSize);
            filter = 1./(ellipse.^factor);       % Construct the filter.
            filter = fftshift(filter);
            phasemod = fftshift(ellipse.^velfactor + lowvel);
        end
      end
      
      dphase = random('norm',meandev,stddev,filterSize,filterSize);  % Random normal distribution 
      dphase = dphase.*phasemod;

      phase     = phase + dphase;
      newfft    = filter .* exp(1i*phase);
      im        = real(ifft2(newfft));          % Invert to obtain image in spatial domain

      colOffset = (size(im,2) - desiredMatrixSize(2))/2;
      rowOffset = (size(im,1) - desiredMatrixSize(1))/2;
      im        = im(rowOffset + (1:desiredMatrixSize(1)), colOffset + (1:desiredMatrixSize(2)));
     
      % normalized to 1.0
      im = im / max(abs(im(:)));
      
      % zero mean
      im = im - mean(im(:));
      
      threeLevels = false;
      if (threeLevels)
          im(find(im <-1)) = -1;
            im(find(im > 1)) = 1;
          % make it a 3-level image (-1, 0, 1)
          im(find(abs(im) < 0.33)) = 0;
          im(find(im < -eps)) = -1;
          im(find(im > eps)) = 1;
      else
          % two levels
          negIndicies = find(im < -eps);
          im = ones(size(im));
          im(negIndicies) = -1;
      end

      % normalize to [0 .. 1];
      im = (im + 1)/2;
      
      % store sequence
      originalIm  = uint8(255*im);
      imageSequence(1, frameIndex,:,:) = originalIm;
      
      % make pixelated versions
      for blockSizeIndex = 1:numel(blockSizeArray)
          blockSize = blockSizeArray(blockSizeIndex);
            pixelatedIm = pixelateImage(im, blockSize);
            pixelatedIm = uint8(255*pixelatedIm);
            [mean(originalIm(:)) mean(pixelatedIm(:))]
            imageSequence(1+blockSizeIndex, frameIndex,:,:) = pixelatedIm;
      end
      
    end % frameIndex
end


function pixelatedIm = pixelateImage(im, blockSize)
    pixelatedIm = zeros(size(im)); 
    for row = 1:blockSize:size(im,1)
        minRow = row;
        maxRow = min([minRow+blockSize-1 size(im,1)]);
        rowIndices = minRow:maxRow;
        for col = 1:blockSize:size(im,2)
            minCol = col;
            maxCol = min([minCol + blockSize-1 size(im,2)]);
            colIndices = minCol:maxCol;
            dataToAverage = im(rowIndices, colIndices);
            pixelatedIm(rowIndices, colIndices) = mean(dataToAverage(:));
        end
    end 
end


function filter = makeEllipseFilter(orientation, ratio, size)
    angle = orientation/180*pi;
    x = [-size/2:size/2-1];
    y = x;
    [X,Y] = meshgrid(x,y);
    XX = X * cos(angle) - Y*sin(angle);
    YY = X * sin(angle) + Y*cos(angle);
    filter = sqrt((XX).^2 + (YY/ratio).^2);
    filter(size/2+1,size/2+1) = 1;      % .. avoid division by zero.
end