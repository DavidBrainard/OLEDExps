function OLDAnalyzePreliminarySamsungOLEDdata

    close all
    clear all
    clear classes
    clc
    
    % Set following flag to true to export all calibration frames as pdfs
    printCalibrationFrames = false;
    
    
    % Data file where all data structs are appended
    % calibrationFileName = '/Users/Shared/Matlab/Toolboxes/BrainardLabToolbox/OOCalibrationToolbox/SamsungOLED_calib.mat';
    
    %  Single target runs expect for last 2
    calibrationFileName = './PreliminaryData/SamsungOLED_calib.mat';
    
    % Double target runs
    calibrationFileName = './PreliminaryData/SamsungOLED_DoubleTargetCalib1.mat';
    
    % create a MAT-file object that supports partial loading and saving.
    matOBJ = matfile(calibrationFileName, 'Writable', false);
    
    % get current variables
    varList = who(matOBJ);
        
    if isempty(varList)
        if (exist(dataSetFilename, 'file'))
            fprintf(2,'No calibration data found in ''%s''.\n', dataSetFilename);
        else
            fprintf(2,'''%s'' does not exist.\n', dataSetFilename);
        end
        calibrationDataSet = [];
        return;        
    end
    
    fprintf('\nFound %d calibration data sets in the saved history.', numel(varList));
    
    % ask the user to select one
    defaultDataSetNo = numel(varList);
    dataSetIndex = input(sprintf('\nSelect a data set (1-%d) [%d]: ', defaultDataSetNo, defaultDataSetNo));
    if isempty(dataSetIndex) || (dataSetIndex < 1) || (dataSetIndex > defaultDataSetNo)
       dataSetIndex = defaultDataSetNo;
    end
      
    % return the selected ground truth data set
    eval(sprintf('calibrationDataSet = matOBJ.%s;',varList{dataSetIndex}));
    
    % Retrieve data
    allCondsData = calibrationDataSet.allCondsData;
    runParams = calibrationDataSet.runParams
    
    stabilizerBorderWidth   = runParams.stabilizerBorderWidth;
    stabilizerGrayLevelNum  = numel(runParams.stabilizerGrays);
    sceneGrayLevelNum       = numel(runParams.sceneGrays);
    biasLevelNum            = numel(runParams.biasGrays);
    biasSizesNum            = size(runParams.biasSizes,1);
    gammaInputValuesNum     = numel(runParams.leftTargetGrays);
    
    % Load CIE 1931 CMFs
    load T_xyz1931
    vLambda1931_originalSampling = squeeze(T_xyz1931(2,:));
    desiredS = [380 1 401];
    
    
    
    fprintf('\n\n'); 
    fprintf('\n%-30s: %s', 'Temporal dithering mode', runParams.temporalDitheringMode);
    
    fprintf('\n%-30s: %d pixels', 'Stabilizer border width', runParams.stabilizerBorderWidth);
   
    fprintf('\n%-30s: ', 'Stabilizer gray level(s)');
    fprintf('%2.2f ', runParams.stabilizerGrays);
    
    fprintf('\n%-30s: ', 'Scene mean gray level(s)');
    fprintf('%2.2f ', runParams.sceneGrays);
    
    fprintf('\n%-30s: ', 'Bias gray level(s)');
    fprintf('%2.2f ', runParams.biasGrays);
    
    fprintf('\n%-30s: ', 'Bias region sizes (x)');
    fprintf('%2.2f ', runParams.biasSizes(:, 1));
    
    fprintf('\n%-30s: ', 'Bias region sizes (y)');
    fprintf('%2.2f ', runParams.biasSizes(:, 2));
    
    fprintf('\n%-30s: ', 'Gamma input values (left)');
    fprintf('%2.3f ', runParams.leftTargetGrays);
    
    fprintf('\n%-30s: ', 'Gamma input values (right)');
    fprintf('%2.3f ', runParams.rightTargetGrays);
    
    fprintf('\n%-30s: (%d,%d)', 'Target position (left)', runParams.leftTarget.x0, runParams.leftTarget.y0);
    fprintf('\n%-30s: (%d,%d)', 'Target position (right)', runParams.rightTarget.x0, runParams.rightTarget.y0);
    
    fprintf('\n%-30s: %dx%d', 'Target size (left)', runParams.leftTarget.width, runParams.leftTarget.height);
    fprintf('\n%-30s: %dx%d', 'Target size (right)', runParams.rightTarget.width, runParams.rightTarget.height);
    
    
    
    fprintf('\n\n');
    
    
    % Preallocate memory for spds
    leftSPD = zeros(stabilizerGrayLevelNum, ...
                    sceneGrayLevelNum, ...
                    biasLevelNum, ...
                    biasSizesNum, ...
                    gammaInputValuesNum, ...
                    desiredS(3));
           
    rightSPD = leftSPD;
    
    
    conditionsNum = numel(allCondsData);
    
    
    totalEnergy = zeros(stabilizerGrayLevelNum, biasSizesNum,gammaInputValuesNum);

    for condIndex = 1:conditionsNum
        
        condIndex
        % get struct for current condition
        conditionData = allCondsData{condIndex};
       
        % get indices
        stabilizerGrayIndex = conditionData.stabilizerGrayIndex;
        sceneGrayIndex = conditionData.sceneGrayIndex;
        biasGrayIndex = conditionData.biasGrayIndex;
        biasSizeIndex = conditionData.biasSizeIndex;
        leftTargetGrayIndex = conditionData.leftTargetGrayIndex;
        rightTargetGrayIndex = conditionData.rightTargetGrayIndex;
        
        
        % analyze the stimFrame to extract various energies
        stimFrame = double(squeeze(conditionData.demoFrame(:,:,1)));
        stimFrame(find(stimFrame<0)) = 0;
        stimFrame(find(stimFrame>1)) = 1;
    
        
        
        [biasRegionIndices, leftTargetRegionIndices, rightTargetRegionIndices, sceneRegionIndices, stabilizerRegionIndices] = ...
            computeComponentIndices(stimFrame, biasSizeIndex, runParams); 
        
        
        biasRegionEnergy(stabilizerGrayIndex, biasSizeIndex, leftTargetGrayIndex)        = sum(stimFrame(biasRegionIndices));
        leftTargetRegionEnergy(stabilizerGrayIndex, biasSizeIndex, leftTargetGrayIndex)  = sum(stimFrame(leftTargetRegionIndices));
        rightTargetRegionEnergy(stabilizerGrayIndex, biasSizeIndex, leftTargetGrayIndex) = sum(stimFrame(rightTargetRegionIndices));
        sceneRegionEnergy(stabilizerGrayIndex, biasSizeIndex, leftTargetGrayIndex)       = sum(stimFrame(sceneRegionIndices));
        stabilizerRegionEnergy(stabilizerGrayIndex, biasSizeIndex, leftTargetGrayIndex)  = sum(stimFrame(stabilizerRegionIndices));
        
        showRegionAnalysisResults = false;
        if (showRegionAnalysisResults)
            
             % just for testing
%             stimFrame(biasRegionIndices) = 0.5;
%             stimFrame(leftTargetRegionIndices) = 0.25;
%             stimFrame(rightTargetRegionIndices) = 0.75;
%             stimFrame(sceneRegionIndices) = 0.9;
%             stimFrame(stabilizerRegionIndices) = 0.2;
        
            [rows, cols] = size(stimFrame);
            figure(99);
            clf;
            gain = rows/4;
            for currentRow = 1:rows
                
                sliceProfile = squeeze(stimFrame(currentRow,:));
                
                subplot(1,2,1);
                imagesc([1:cols], [1:rows], stimFrame);
                hold on;
                stairs([1:cols], currentRow - sliceProfile*gain, 'r-');
                hold off;
                colormap(gray);
                set(gca, 'CLim', [0 1]);
                axis 'ij'
                
                drawnow;
                (0.2);
            end

            subplot(1,2,2);
            stairs([1:cols], double(stimFrame(rows/2,:)), 'r-');
            set(gca, 'YLim', [0 1]);
        end % showRegionAnalysisResults
    
        
        if (printCalibrationFrames)
            h0 = figure(99);
            set(h0, 'Position', [100 100 754 453]);
            imshow(stimFrame);
            hold on;
            plot([1 size(stimFrame,2) size(stimFrame,2) 1 1], ...
                 [1 1 size(stimFrame,1) size(stimFrame,1)  1], 'k-');
            hold off;
            colormap(gray(256));
            set(gca, 'CLim', [0 1]);
            axis 'image'
            drawnow;

            % Print frame as pdf
            set(h0,'PaperOrientation','landscape');
            set(h0,'PaperUnits','normalized');
            set(h0,'PaperPosition', [0 0 1 1]);
            print(gcf, '-dpdf', sprintf('Cond_%d.pdf', condIndex));
        end % printCalibrationFrames
        
    
        if (condIndex == 1)
            nativeS = conditionData.leftS;
            vLambda = 683*SplineCmf(S_xyz1931, vLambda1931_originalSampling, desiredS);
            wave = SToWls(desiredS);
        end
        
        
        totalEnergy(stabilizerGrayIndex, biasSizeIndex, leftTargetGrayIndex) = sum(sum(stimFrame));
        leftGammaIn(stabilizerGrayIndex, biasSizeIndex, leftTargetGrayIndex) = runParams.leftTargetGrays(leftTargetGrayIndex);
        rightGammaIn(stabilizerGrayIndex, biasSizeIndex, rightTargetGrayIndex) = runParams.rightTargetGrays(rightTargetGrayIndex);
        
        
        % get SPD data 
        spd = conditionData.leftSPD;
        
        % interpolate to desiredS
        spd = SplineSpd(nativeS, spd', desiredS);
        
        leftSPD(stabilizerGrayIndex, ...
                sceneGrayIndex, ...
                biasGrayIndex, ...
                biasSizeIndex, ...
                leftTargetGrayIndex, ...
                :) = spd;
        

        leftGammaOut(stabilizerGrayIndex, biasSizeIndex, leftTargetGrayIndex) = sum(spd'.*vLambda);
            
        if ~isempty(conditionData.rightSPD) 
            % get SPD data 
            spd = conditionData.rightSPD;
            % interpolate to desiredS
            spd = SplineSpd(nativeS, spd', desiredS);
        
            rightSPD(stabilizerGrayIndex, ...
                    sceneGrayIndex, ...
                    biasGrayIndex, ...
                    biasSizeIndex, ...
                    rightTargetGrayIndex, ...
                    :) = spd;
            rightGammaOut(stabilizerGrayIndex, biasSizeIndex, leftTargetGrayIndex) = sum(spd'.*vLambda);
        else % ~isempty(conditionData.rightSPD)
            rightSPD = [];
            rightGammaOut = [];
        end
    end % cond Index
    
    
    size(leftGammaOut)
    size(rightGammaOut)
    
    
    %totalEnergy(stabilizerGrayIndex, biasSizeIndex, leftTargetGrayIndex) 
    %leftGammaIn(stabilizerGrayIndex, biasSizeIndex, leftTargetGrayIndex) 
    %leftGammaOut(stabilizerGrayIndex, biasSizeIndex, leftTargetGrayIndex) 
    %rightGammaIn(stabilizerGrayIndex, biasSizeIndex, rightTargetGrayIndex)
    %rightGammaOut(stabilizerGrayIndex, biasSizeIndex, leftTargetGrayIndex)  
    
    h99 = figure(99);
    clf;
    subplot('Position', [0.04 0.04 0.9 0.9]);
    referenceGammaCurveLeft = squeeze(leftGammaOut(1, 1,:));
    referenceGammaCurveRight = squeeze(rightGammaOut(1, 1,:));
    
    hold on
    for stabilizerGrayIndex = 1:stabilizerGrayLevelNum
        for biasSizeIndex = 1: biasSizesNum
            gammaCurveLeft = squeeze(leftGammaOut(stabilizerGrayIndex, biasSizeIndex,:));
            gammaCurveRight = squeeze(rightGammaOut(stabilizerGrayIndex, biasSizeIndex,:));
             
            ratiosLeft(stabilizerGrayIndex, biasSizeIndex,:) = gammaCurveLeft./referenceGammaCurveLeft;
            ratiosRight(stabilizerGrayIndex, biasSizeIndex,:) = gammaCurveRight./referenceGammaCurveRight;
            curveRatiosLeft(stabilizerGrayIndex, biasSizeIndex) = gammaCurveLeft \ referenceGammaCurveLeft;
            curveRatiosRight(stabilizerGrayIndex, biasSizeIndex) = gammaCurveRight \ referenceGammaCurveRight;
            
            
            plot(squeeze(totalEnergy(stabilizerGrayIndex, biasSizeIndex,:)), gammaCurveLeft, 'ro-');
            plot(squeeze(totalEnergy(stabilizerGrayIndex, biasSizeIndex,:)), gammaCurveRight, 'bo-');
            
        end
    end
    hold off;
    
    xlabel('Stimulus energy');
    set(gca, 'XLim', [min(totalEnergy(:)) max(totalEnergy(:))]);
    
    
    h100 = figure(100);
    clf;
    subplot('Position', [0.04 0.04 0.9 0.9]);
    
    hold on
    for stabilizerGrayIndex = 1:stabilizerGrayLevelNum
        for biasSizeIndex = 1: biasSizesNum
            
            
            plot(squeeze(totalEnergy(stabilizerGrayIndex, biasSizeIndex,:)), repmat(squeeze(curveRatiosLeft(stabilizerGrayIndex, biasSizeIndex)), [1 size(totalEnergy,3)]), 'r-');
            plot(squeeze(totalEnergy(stabilizerGrayIndex, biasSizeIndex,:)), repmat(squeeze(curveRatiosRight(stabilizerGrayIndex, biasSizeIndex)), [1 size(totalEnergy,3)]), 'b-');
            
            indices = find(~isnan(squeeze(ratiosLeft(stabilizerGrayIndex, biasSizeIndex,:))));
            plot(squeeze(totalEnergy(stabilizerGrayIndex, biasSizeIndex,indices)), squeeze(ratiosLeft(stabilizerGrayIndex, biasSizeIndex, indices)), 'r-');
            indices = find(~isnan(squeeze(ratiosRight(stabilizerGrayIndex, biasSizeIndex,:))));
            plot(squeeze(totalEnergy(stabilizerGrayIndex, biasSizeIndex,indices)), squeeze(ratiosRight(stabilizerGrayIndex, biasSizeIndex, indices)), 'b-');
            
           
        end
    end
    hold off;
     
    xlabel('Stimulus energy');
    ylabel('Output Scaling');
    set(gca, 'XLim', [min(totalEnergy(:)) max(totalEnergy(:))]);
    
    

    
    
    % plot data 
    
    stabilizerGrayIndex = 1;
    sceneGrayIndex = 1;
    biasGrayIndex = 1;
    biasSizeIndex = 1;
    
    stabilizerGray = runParams.stabilizerGrays(stabilizerGrayIndex);
    sceneGray      = runParams.sceneGrays(sceneGrayIndex);
    biasGray       = runParams.biasGrays(biasGrayIndex);
    
    
    gammaOutputLeft = zeros(stabilizerGrayLevelNum, biasSizesNum, gammaInputValuesNum);
    gammaOutputRight = zeros(stabilizerGrayLevelNum, biasSizesNum, gammaInputValuesNum);
    vLambda  = repmat(vLambda, [gammaInputValuesNum 1]);
    
    for stabilizerGrayIndex = 1:stabilizerGrayLevelNum
        for biasSizeIndex = 1: biasSizesNum
            spd = squeeze(leftSPD(stabilizerGrayIndex, ...
                        sceneGrayIndex, ...
                        biasGrayIndex, ...
                        biasSizeIndex, ...
                        1:gammaInputValuesNum, ...
                        :));
            luminance = sum(spd.*vLambda,2);
            gammaOutputLeft(stabilizerGrayIndex, biasSizeIndex, :) = luminance; 
        end
    end
    
    
    
    if ~isempty(rightSPD)
        for stabilizerGrayIndex = 1:stabilizerGrayLevelNum
            for biasSizeIndex = 1: biasSizesNum
                spd = squeeze(rightSPD(stabilizerGrayIndex, ...
                            sceneGrayIndex, ...
                            biasGrayIndex, ...
                            biasSizeIndex, ...
                            1:gammaInputValuesNum, ...
                            :));
                luminance = sum(spd.*vLambda,2);
                gammaOutputRight(stabilizerGrayIndex, biasSizeIndex, :) = luminance; 
            end
        end
    end
    
    
    
        
        
    
    
    
    gammaInputLeft  = runParams.leftTargetGrays;
    maxGammaOutputLeft = max(gammaOutputLeft(:));
    
    if isempty(rightSPD)
        gammaInputRight = [];
    else
        gammaInputRight  = runParams.rightTargetGrays;
        maxGammaInputRight = max(gammaInputRight(:));
        minGammaInputRight = min(gammaInputRight(:));
    end
    
    
    
    
    
    h1 = figure(1);
    figXo = 2560;   figYo = 360;
    figWidth = 700; figHeight = 860;
    set(h1, 'Position', [figXo figYo figWidth figHeight]);
    clf;
    
    lineColors = lines(stabilizerGrayLevelNum*biasSizesNum);
    
    % Subplot panel sizes and margins
    width = 0.85/(biasSizesNum+1);
    height = 0.7/(stabilizerGrayLevelNum+1);
    marginX = 0.02;
    marginY = 0.05;
    
	
    referenceBiasSizeIndex = 1;
    referenceBiasSizeX  = runParams.biasSizes(referenceBiasSizeIndex,1);
    referenceBiasSizeY  = runParams.biasSizes(referenceBiasSizeIndex,2);
        
    % First scan (all condition curves plus scaled curves for left to right conditions)
    for stabilizerGrayIndex = 1:stabilizerGrayLevelNum
        
        stabilizerGray      = runParams.stabilizerGrays(stabilizerGrayIndex);
        referenceGammaCurve = squeeze(gammaOutputLeft(stabilizerGrayIndex, referenceBiasSizeIndex,:));
        
        legendMatrix = {};
        for biasSizeIndex = 1: biasSizesNum
            
            biasSizeX = runParams.biasSizes(biasSizeIndex, 1);
            biasSizeY = runParams.biasSizes(biasSizeIndex, 2);
        
            gammaCurveLeft   = squeeze(gammaOutputLeft(stabilizerGrayIndex, biasSizeIndex, :));
            scalingFactor    = gammaCurveLeft \ referenceGammaCurve;
            scaledGammaCurve = gammaCurveLeft * scalingFactor;
            
            left = 3*marginX + (biasSizeIndex-1)*(width+marginX);
            bottom = 1-stabilizerGrayIndex*(height+marginY);
            subplot('Position', [left bottom width height]);   
            
            condIndex = (stabilizerGrayIndex-1)* biasSizesNum + biasSizeIndex;
            lineColor = lineColors(condIndex,:);
            
            plot(gammaInputLeft, gammaCurveLeft, 'ks-', 'LineWidth', 3.0, 'MarkerSize', 8, 'MarkerFaceColor', [0.8 0.8 0.8], 'Color', lineColor);
            
            if (~isempty(gammaInputRight))
                gammaCurveRight = squeeze(gammaOutputRight(stabilizerGrayIndex, biasSizeIndex, :));
                if (maxGammaInputRight == minGammaInputRight)
                    gammaInput = gammaInputLeft;
                else
                    gammaInput = gammaInputRight;
                end
                
                hold on;
                plot(gammaInput, gammaCurveRight, 'k.-', 'LineWidth', 3.0, 'MarkerSize', 8);
                if (maxGammaInputRight == minGammaInputRight)
                    plot([minGammaInputRight, maxGammaInputRight], [0 mean(gammaCurveRight)], 'k-', 'LineWidth', 3.0);
                    plot(minGammaInputRight, 0, 'kv', 'LineWidth', 2.0, 'MarkerFaceColor', [0.6 0.6 0.6]);
                end
                hold off;
                legend_handle = legend({'Left', 'Right'}, 'FontName', 'Helvetica', 'FontSize', 6, 'Location', 'NorthWest');
                set(legend_handle, 'Box', 'off')
            end
            
            set(gca, 'FontName', 'Helvetica', 'FontSize', 8);
            grid on;
            box on
            
            set(gca, 'YLim', [0 maxGammaOutputLeft]);
            set(gca, 'XTick', [0:0.2:1.0], 'YTick', [0:100:1000]);
            
            if (biasSizeIndex == 1)
               ylabel('luminance (cd/m2)', 'FontName', 'Helvetica', 'FontSize', 10, 'FontWeight', 'bold');
               set(gca, 'YTickLabel', [0:100:1000]);
            else
               ylabel(''); 
               set(gca, 'YTickLabel', []);
            end
            xlabel('');
            
            title(sprintf('Stabilizer gray = %2.2f; \nBias WxH = %2.0fx%2.0f pxls.', stabilizerGray, biasSizeX, biasSizeY), 'FontName', 'Helvetica', 'FontSize', 8);
            
            % The scaled gamma curves for the CurrentStabilizerGray
            left = 3*marginX + biasSizesNum*(width+marginX);
            subplot('Position', [left bottom width height]);
            
            hold on;
            plot(gammaInputLeft, scaledGammaCurve, 'k-', 'LineWidth', 3.0, 'MarkerSize', 8, 'MarkerFaceColor', [0.8 0.8 0.8], 'Color', lineColors(condIndex,:));
            legendMatrix{biasSizeIndex} = sprintf('BiasWxH: %2.0fx%2.0f (scale: %2.2f)', biasSizeX, biasSizeY, 1.0/scalingFactor);     
        end
        
        xlabel('');
        ylabel('');
        set(gca, 'YTickLabel', []);
        set(gca, 'YLim', [0 maxGammaOutputLeft]);
        set(gca, 'XTick', [0:0.2:1.0], 'YTick', [0:100:1000]);
        set(gca, 'FontName', 'Helvetica', 'FontSize', 8);
        grid on;
        box on
            
        % legend and title
        legend_handle = legend(legendMatrix, 'FontName', 'Helvetica', 'FontSize', 6, 'Location', 'Best');
        set(legend_handle, 'Box', 'off')
        title(sprintf('Scaled gammas w/r to:\nBiasWxH = %2.2f x %2.2f pxls', referenceBiasSizeX, referenceBiasSizeY), 'FontName', 'Helvetica', 'FontSize', 8, 'BackgroundColor',[.99 .99 .48], 'EdgeColor', [0 0 0]);
    end
    
   
    % Second scan  (scaled curves for top-to-bottom conditions)
    referenceStabilizerGrayIndex = 1;
    for biasSizeIndex = 1: biasSizesNum
        referenceGammaCurve = squeeze(gammaOutputLeft(referenceStabilizerGrayIndex, biasSizeIndex,:));
        referenceStabilizerGray  = runParams.stabilizerGrays(referenceStabilizerGrayIndex);
        
        legendMatrix = {};
        for stabilizerGrayIndex = 1:stabilizerGrayLevelNum
            stabilizerGray   = runParams.stabilizerGrays(stabilizerGrayIndex);
            gammaCurve       = squeeze(gammaOutputLeft(stabilizerGrayIndex, biasSizeIndex, :));
            scalingFactor    = gammaCurve \ referenceGammaCurve;
            scaledGammaCurve = gammaCurve * scalingFactor;
            
            left = 3*marginX + (biasSizeIndex-1)*(width+marginX);
            bottom = 1-(stabilizerGrayLevelNum+1)*(height+marginY);
            subplot('Position', [left bottom width height]);   
            
            condIndex = (stabilizerGrayIndex-1)* biasSizesNum + biasSizeIndex;
            lineColor = lineColors(condIndex,:);
            
            hold on;
            plot(gammaInputLeft, scaledGammaCurve, 'k-', 'LineWidth', 3.0, 'MarkerSize', 8, 'MarkerFaceColor', [0.8 0.8 0.8], 'Color', lineColors(condIndex,:));
            legendMatrix{stabilizerGrayIndex} = sprintf('Stabil. gray = %2.2f (scale: %2.2f)', stabilizerGray, 1.0/scalingFactor);
        end
        
        xlabel('settings value', 'FontName', 'Helvetica', 'FontSize', 10, 'FontWeight', 'bold');
        if (biasSizeIndex == 1)
               ylabel('luminance (cd/m2)', 'FontName', 'Helvetica', 'FontSize', 10, 'FontWeight', 'bold');
               set(gca, 'YTickLabel', [0:100:1000]);
        else
               ylabel(''); 
               set(gca, 'YTickLabel', []);
        end
            
        set(gca, 'YLim', [0 maxGammaOutputLeft]);
        set(gca, 'XTick', [0:0.2:1.0], 'YTick', [0:100:1000]);
        set(gca, 'FontName', 'Helvetica', 'FontSize', 8);
        grid on;
        box on
        
        % legend and title
        legend_handle = legend(legendMatrix, 'FontName', 'Helvetica', 'FontSize', 6, 'Location', 'Best');
        set(legend_handle, 'Box', 'off')
        title(sprintf('Scaled gammas w/r to:\nStabilizerGray = %2.2f', referenceStabilizerGray), 'FontName', 'Helvetica', 'FontSize', 8, 'BackgroundColor',[.99 .99 .48], 'EdgeColor', [0 0 0]);
    end % biasSizeIndex
    

    % Third scan (all scaled curves)
    referenceGammaCurve = squeeze(gammaOutputLeft(referenceStabilizerGrayIndex, referenceBiasSizeIndex,:));
    referenceStabilizerGray   = runParams.stabilizerGrays(referenceStabilizerGrayIndex);
            
    condIndex = 0;
    legendMatrix = {};
    
    left = 3*marginX + biasSizesNum*(width+marginX);
    bottom = 1-(stabilizerGrayLevelNum+1)*(height+marginY);
    subplot('Position', [left bottom width height]);   
    hold on;
    
    for stabilizerGrayIndex = 1:stabilizerGrayLevelNum
        for biasSizeIndex = 1: biasSizesNum
            gammaCurve       = squeeze(gammaOutputLeft(stabilizerGrayIndex, biasSizeIndex, :));
            scalingFactor    = gammaCurve \ referenceGammaCurve;
            scaledGammaCurve = gammaCurve * scalingFactor;
            scalingFactorMatrix(stabilizerGrayIndex, biasSizeIndex) = 1.0/scalingFactor;
            
            condIndex = (stabilizerGrayIndex-1)* biasSizesNum + biasSizeIndex;
            lineColor = lineColors(condIndex,:);
            plot(gammaInputLeft, scaledGammaCurve, 'k-', 'LineWidth', 3.0, 'MarkerSize', 8, 'MarkerFaceColor', [0.8 0.8 0.8], 'Color', lineColors(condIndex,:));
            legendMatrix{condIndex} = sprintf('scale: %2.2f', 1.0/scalingFactor);
        end       
    end % biasSizeIndex
    
    xlabel('settings value', 'FontName', 'Helvetica', 'FontSize', 10, 'FontWeight', 'bold');
    ylabel('');
    set(gca, 'YTickLabel', []);
    set(gca, 'YLim', [0 maxGammaOutputLeft]);
    set(gca, 'XTick', [0:0.2:1.0], 'YTick', [0:100:1000]);
    set(gca, 'FontName', 'Helvetica', 'FontSize', 8);
    grid on;
    box on
    
    % legend and title
    legend_handle = legend(legendMatrix, 'FontName', 'Helvetica', 'FontSize', 6, 'Location', 'NorthWest');
    set(legend_handle, 'Box', 'off')
    title(sprintf('Scaled gammas w/r to:\nStab.Gray=%2.2f, BiasWxH=%2.0fx%2.0f pxls.', referenceStabilizerGray, referenceBiasSizeX, referenceBiasSizeY), 'FontName', 'Helvetica', 'FontSize', 8, 'BackgroundColor',[.99 .99 .48], 'EdgeColor', [0 0 0]);

        
    
    % Print figure
    set(h1,'PaperOrientation','Portrait');
    set(h1,'PaperUnits','normalized');
    set(h1,'PaperPosition', [0 0 1 1]);
    print(gcf, '-dpdf', '-r600', 'Fig1.pdf');
    
    
    figure(2);
    imagesc(scalingFactorMatrix);
    colormap(gray);
    colorbar
    set(gca, 'CLim', [0.5 1]);
    drawnow;
    
    
    figure(3);
    clf;

    
    for stabilizerGrayIndex = 1:stabilizerGrayLevelNum
        for biasSizeIndex = 1: biasSizesNum 
            
            luminanceLeft  = squeeze(gammaOutputLeft(stabilizerGrayIndex, biasSizeIndex, :));
            luminanceRight  = squeeze(gammaOutputRight(stabilizerGrayIndex, biasSizeIndex, :));
            
            
            biasEnergy = squeeze(biasRegionEnergy(stabilizerGrayIndex, biasSizeIndex, :));
            leftTargetEnergy = squeeze(leftTargetRegionEnergy(stabilizerGrayIndex, biasSizeIndex, :));
            rightTargetEnergy = squeeze(rightTargetRegionEnergy(stabilizerGrayIndex, biasSizeIndex, :));
            sceneEnergy = squeeze(sceneRegionEnergy(stabilizerGrayIndex, biasSizeIndex, :));
            stabilizerEnergy = squeeze(stabilizerRegionEnergy(stabilizerGrayIndex, biasSizeIndex, :));
            
            % mean over gamma input curves, so these should be constant
            biasEnergy = mean(biasEnergy);
            sceneEnergy = mean(sceneEnergy);
            stabilizerEnergy = mean(stabilizerEnergy);
            
            scaleFactor = scalingFactorMatrix(stabilizerGrayIndex, biasSizeIndex);
            subplot(1,4,1);
            hold on;
            plot(stabilizerEnergy+biasEnergy+sceneEnergy, scaleFactor, 'ks');
            
            subplot(1,4,2);
            hold on;
            plot(stabilizerEnergy, scaleFactor, 'ks');
            
            subplot(1,4,3);
            hold on;
            plot(biasEnergy, scaleFactor, 'ks');
            
            subplot(1,4,4);
            hold on;
            plot(sceneEnergy, scaleFactor, 'ks');
            
        end
    end
    
    subplot(1,4,1);
    xlabel('total energy');
    
    subplot(1,4,2);
    xlabel('stabilizer energy');
    
    subplot(1,4,3);
    xlabel('bias energy');
    
    subplot(1,4,4);
    xlabel('scene energy');
    
    
    
    
    
end


function [biasRegionExcludingTargetRegionIndices, leftTargetRegionIndices, rightTargetRegionIndices, sceneRegionExcludingBiasAndTargetRegionIndices, stabilizerRegionIndices] = ...
            computeComponentIndices(stimFrame, biasSizeIndex, runParams)
        % indices of various regions
        biasRegionXindices  = round(0.5*(runParams.leftTarget.x0 - runParams.biasSizes(biasSizeIndex, 1)/2 + [1:runParams.biasSizes(biasSizeIndex, 1)]));
        biasRegionYindices  = round(0.5*(runParams.leftTarget.y0 - runParams.biasSizes(biasSizeIndex, 2)/2 + [1:runParams.biasSizes(biasSizeIndex, 2)]));
        leftTargetXindices  = round(0.5*(runParams.leftTarget.x0 - runParams.leftTarget.width/2 + [1:runParams.leftTarget.width]));
        leftTargetYindices  = round(0.5*(runParams.leftTarget.y0 - runParams.leftTarget.height/2 + [1:runParams.leftTarget.height]));
        rightTargetXindices = round(0.5*(runParams.rightTarget.x0 - runParams.rightTarget.width/2 + [1:runParams.rightTarget.width]));
        rightTargetYindices = round(0.5*(runParams.rightTarget.y0 - runParams.rightTarget.height/2 + [1:runParams.rightTarget.height]));
        sceneRegionXindices = round(0.5*([runParams.stabilizerBorderWidth + [1: size(stimFrame,2)*2-runParams.stabilizerBorderWidth*2]]));
        sceneRegionYindices = round(0.5*([runParams.stabilizerBorderWidth + [1: size(stimFrame,1)*2-runParams.stabilizerBorderWidth*2]]));
        frameXindices    = [1:size(stimFrame,2)];
        frameYindices    = [1:size(stimFrame,1)];
        
        [X,Y] = meshgrid(frameYindices, frameXindices);
        frameIndices = sub2ind(size(stimFrame), X,Y);
  
        
        [X,Y] = meshgrid(biasRegionYindices, biasRegionXindices);
        biasRegionIndices = sub2ind(size(stimFrame), X,Y);
        
        [X,Y] = meshgrid(leftTargetYindices, leftTargetXindices);
        leftTargetRegionIndices = sub2ind(size(stimFrame), X,Y);
        
        [X,Y] = meshgrid(rightTargetYindices, rightTargetXindices);
        rightTargetRegionIndices = sub2ind(size(stimFrame), X,Y);
        
        [X,Y] = meshgrid(sceneRegionYindices, sceneRegionXindices);
        sceneRegionIndices = sub2ind(size(stimFrame), X,Y);
        
        stabilizerRegionIndices = setdiff(frameIndices, sceneRegionIndices);
        
        biasRegionExcludingTargetRegionIndices = setdiff(biasRegionIndices, leftTargetRegionIndices);
        biasRegionExcludingTargetRegionIndices = setdiff(biasRegionExcludingTargetRegionIndices, rightTargetRegionIndices);
        
        sceneRegionExcludingBiasAndTargetRegionIndices = setdiff(sceneRegionIndices, biasRegionIndices);
        sceneRegionExcludingBiasAndTargetRegionIndices = setdiff(sceneRegionExcludingBiasAndTargetRegionIndices, leftTargetRegionIndices);
        sceneRegionExcludingBiasAndTargetRegionIndices = setdiff(sceneRegionExcludingBiasAndTargetRegionIndices, rightTargetRegionIndices);
        
        biasRegionExcludingTargetRegionIndices = biasRegionExcludingTargetRegionIndices(:);
        leftTargetRegionIndices = leftTargetRegionIndices(:);
        rightTargetRegionIndices = rightTargetRegionIndices(:);
        sceneRegionExcludingBiasAndTargetRegionIndices = sceneRegionExcludingBiasAndTargetRegionIndices(:);
        stabilizerRegionIndices = stabilizerRegionIndices(:);
        
        
end
