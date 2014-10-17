function playSequence

    load('PixelOLEDprobes.mat');  %loads 'stimParams', 'stimuli'
     
    % start a video writer
    writerObj = VideoWriter('SVMtest.mp4', 'MPEG-4'); 
    writerObj.FrameRate = 30;
    writerObj.Quality = 100;
    open(writerObj);
    
    for exponentOfOneOverFIndex = 1:numel(stimParams.exponentOfOneOverFArray)
        for oriBiasIndex = 1:numel(stimParams.oriBiasArray)
            sequence = stimuli{exponentOfOneOverFIndex, oriBiasIndex}.imageSequence;
            for frameIndex = 1:size(sequence, 2)
                image = squeeze(sequence(1, frameIndex,:,:));
                for pixelSizeIndex = 1:numel(stimParams.blockSizeArray)
                    pixelatedImage(pixelSizeIndex,:,:) = squeeze(sequence(1+pixelSizeIndex, frameIndex,:,:));
                end
                [mean(image(:)) mean(squeeze(pixelatedImage(1,:,:))) mean(squeeze(pixelatedImage(2,:,:))) mean(squeeze(pixelatedImage(3,:,:)))]
                drawFrame(image, pixelatedImage, stimParams.blockSizeArray, writerObj);
            end
        end
    end
    
    
    % now play the inverted polarity
    
     
    % close video writer
    close(writerObj); 
    
end


function drawFrame(image, pixelatedImage, blockSizeArray, writerObj)

    h = figure(1);
    set(h, 'Position', [500 500 931 554]);
    clf;
    colormap(gray(256));
    subplot('Position', [0.03 0.52 0.45 0.45]);
    imagesc(image);
    set(gca, 'CLim', [0 255]);
    hold on;
    xLeft = 1200-10+8;
    yLeft = 650-5+1;
    plot(xLeft, yLeft, 'rs', 'MarkerSize', 30, 'MarkerFaceColor', 'r');
    
    xRight = 1700-10+8;
    yRight = 950-5+1;
    plot(xRight, yRight, 'rs', 'MarkerSize', 30, 'MarkerFaceColor', 'g');
    
    colorbar
    hold off;
    axis 'image'
    axis 'tight'
    set(gca, 'XTick', []);
      
    for blockSizeIndex = 1:numel(blockSizeArray)
        if (blockSizeIndex == 1)
            subplot('Position', [0.52 0.52 0.45 0.45]);
        elseif (blockSizeIndex == 2)
            subplot('Position', [0.03 0.04 0.45 0.45]);
        elseif (blockSizeIndex == 3)
            subplot('Position', [0.52 0.04 0.45 0.45]);
        end
            
        imagesc(squeeze(pixelatedImage(blockSizeIndex,:,:)));
        set(gca, 'CLim', [0 255]);
        hold on;
      
        rowsNum = 1080/blockSizeArray(blockSizeIndex);
        colsNum = 1920/blockSizeArray(blockSizeIndex);
        
        for row = 1:rowsNum
              plot([1 1920], row*blockSizeArray(blockSizeIndex)*[1 1], 'k-');
        end
        for col = 1:colsNum
              plot(col*blockSizeArray(blockSizeIndex)*[1 1], [1 1080],'k-');
        end

        xLeft = 1200-10+8;
        yLeft = 650-5+1;
        plot(xLeft, yLeft, 'rs', 'MarkerSize', 30, 'MarkerFaceColor', 'r');

        xRight = 1700-10+8;
        yRight = 950-5+1;
        plot(xRight, yRight, 'rs', 'MarkerSize', 30, 'MarkerFaceColor', 'g');
    
        colorbar
        hold off;
        axis 'image'
        axis 'tight'
    end
    
    drawnow;
      
    frame = getframe(h);
    writeVideo(writerObj,frame);
end
