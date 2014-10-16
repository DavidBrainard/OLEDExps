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
                pixelatedImage = squeeze(sequence(2, frameIndex,:,:));
                [mean(image(:)) mean(pixelatedImage(:))]
                drawFrame(image, pixelatedImage, stimParams.blockSize, writerObj);
            end
        end
    end
    
    
    % now the inverted polarity
    for exponentOfOneOverFIndex = 1:numel(stimParams.exponentOfOneOverFArray)
        for oriBiasIndex = 1:numel(stimParams.oriBiasArray)
            sequence = stimuli{exponentOfOneOverFIndex, oriBiasIndex}.imageSequence;
            for frameIndex = 1:size(sequence, 2)
                image = 255-squeeze(sequence(1, frameIndex,:,:));
                pixelatedImage = 255-squeeze(sequence(2, frameIndex,:,:));
                -[mean(image(:)) mean(pixelatedImage(:))]
                drawFrame(image, pixelatedImage, stimParams.blockSize, writerObj);
            end
        end
    end
    
     
    % close video writer
    close(writerObj); 
    
end


function drawFrame(image, pixelatedImage, blockSize, writerObj)

    h = figure(1);
    set(h, 'Position', [500 500 790 880]);
    clf;
    colormap(gray(256));
    subplot('Position', [0.03 0.52 0.99 0.45]);
    imagesc(image);
    set(gca, 'CLim', [0 255]);
    hold on;
    for row = 1:11
        plot([1 1920], row*blockSize*[1 1], 'k-');
    end
    for col = 1:20
        plot(col*blockSize*[1 1], [1 1080],'k-');
    end
    x = 1400-10+2;
    y = 820-5+1;
    plot(x,y, 'rs', 'MarkerSize', 30, 'MarkerFaceColor', 'r');
    colorbar
    hold off;
      
    axis 'image'
    axis 'tight'
    set(gca, 'XTick', []);
      
    subplot('Position', [0.03 0.04 0.99 0.45]);
    imagesc(pixelatedImage);
    set(gca, 'CLim', [0 255]);
    hold on;
      
    for row = 1:11
          plot([1 1920], row*blockSize*[1 1], 'k-');
    end
    for col = 1:20
          plot(col*blockSize*[1 1], [1 1080],'k-');
    end

    plot(x,y, 'rs', 'MarkerSize', 30, 'MarkerFaceColor', 'r');
    colorbar
    hold off;
    axis 'image'
    axis 'tight'
    drawnow;
      
    frame = getframe(h);
    writeVideo(writerObj,frame);
end
