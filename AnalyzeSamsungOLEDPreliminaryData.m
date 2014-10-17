function AnalyzePreliminarySamsungOLEDdata

    close all
    clear all
    clear classes
    clc
    
    % Set following flag to true to export all calibration frames as pdfs
    printCalibrationFrames = false;
      
    %  Single target runs expect for last 2 runs
    calibrationFileName = './PreliminaryData/SamsungOLED_calib.mat';
    
    % Double target runs
    %calibrationFileName = './PreliminaryData/SamsungOLED_DoubleTargetCalib1.mat';
    
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
        
        totalEnergy(stabilizerGrayIndex, biasSizeIndex, leftTargetGrayIndex) = sum(sum(stimFrame))/numel(stimFrame);
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
    
    
    
    
    
    h2 = figure(2);
    figXo = 100;   figYo = 100;
    figWidth = 1200; figHeight = 800;
    set(h2, 'Position', [figXo figYo figWidth figHeight]);
    clf;
    
    % The actual gamma curves
    subplot('Position', [0.04 0.74 0.95 0.24]);
    referenceGammaCurveLeft = squeeze(leftGammaOut(1, 1,:));
    referenceGammaCurveRight = referenceGammaCurveLeft; % squeeze(rightGammaOut(1, 1,:));
    
    hold on
    
    lineColors = jet(stabilizerGrayLevelNum*biasSizesNum);
    
    
    for stabilizerGrayIndex = 1:stabilizerGrayLevelNum    
        for biasSizeIndex = 1: biasSizesNum
            
            gammaCurveLeft = squeeze(leftGammaOut(stabilizerGrayIndex, biasSizeIndex,:));
             
            condIndex = (stabilizerGrayIndex-1)* biasSizesNum + biasSizeIndex;
            lineColor = lineColors(condIndex,:);
            plot(squeeze(totalEnergy(stabilizerGrayIndex, biasSizeIndex,:)), gammaCurveLeft, 'o-', 'MarkerSize', 8, 'MarkerFaceColor', lineColor, 'Color', lineColor*0.5, 'LineWidth', 1.0);
            
            if (~isempty(rightGammaOut))
                gammaCurveRight = squeeze(rightGammaOut(stabilizerGrayIndex, biasSizeIndex,:));
                plot(squeeze(totalEnergy(stabilizerGrayIndex, biasSizeIndex,:)), gammaCurveRight, 'ks-', 'MarkerSize', 8, 'MarkerFaceColor', [1 1 1]);
            end
        end
    end
    
    
    energyLims = [0.15 0.25];
    energyMargin = max(totalEnergy(:))*0.01;
    energyLims = [min(totalEnergy(:))-energyMargin max(totalEnergy(:))+energyMargin];
    
    hold off;
    
    set(gca, 'XLim', energyLims, 'YTick', [0:100:1000], 'XTickLabel', []);
    ylabel('luminance (cd/m2)', 'FontName', 'Helvetica', 'FontSize', 10, 'FontWeight', 'bold');
    set(gca, 'FontName', 'Helvetica', 'FontSize', 8, 'Color', [1 1 1]);
    grid on;
   % legend_handle = legend({'left', 'right'}, 'Location', 'NorthEast');
   % set(legend_handle, 'Box', 'on', 'FontName', 'Helvetica', 'FontSize', 8);
    
    box on
    
    
    % The ratios. 
    % Exclude lowest 2 points
    minGammaPointIndexForInclusion = 2;
    indices = minGammaPointIndexForInclusion:gammaInputValuesNum;
    
    subplot('Position', [0.04 0.05 0.95 0.68]);
    hold on
    cond = 0;
    legendMatrix = {};
    for stabilizerGrayIndex = 1:stabilizerGrayLevelNum
        stabilizerGray = runParams.stabilizerGrays(stabilizerGrayIndex);
        for biasSizeIndex = 1: biasSizesNum
            
            condIndex = (stabilizerGrayIndex-1)* biasSizesNum + biasSizeIndex;
            lineColor = lineColors(condIndex,:);
            
            gammaCurveLeft = squeeze(leftGammaOut(stabilizerGrayIndex, biasSizeIndex,:));
            
             % Individual points ratio
            ratiosLeft(stabilizerGrayIndex, biasSizeIndex,:) = gammaCurveLeft./referenceGammaCurveLeft;
            
            if (~isempty(rightGammaOut))
                gammaCurveRight = squeeze(rightGammaOut(stabilizerGrayIndex, biasSizeIndex,:)); 
                ratiosRight(stabilizerGrayIndex, biasSizeIndex,:) = gammaCurveRight./referenceGammaCurveRight;
            end
            
           

             % individual gamma point ratios (left target)
            plot(squeeze(totalEnergy(stabilizerGrayIndex, biasSizeIndex,indices)), squeeze(ratiosLeft(stabilizerGrayIndex, biasSizeIndex, indices)), 'o-', 'MarkerSize', 8, 'MarkerFaceColor', lineColor, 'Color', lineColor*0.5, 'LineWidth', 1); 

            biasSizeX = runParams.biasSizes(biasSizeIndex, 1);
            biasSizeY = runParams.biasSizes(biasSizeIndex, 2);
            cond = cond + 1;
            legendMatrix{cond} = sprintf('BiasWxH: %2.0fx%2.0f, Stabilizer Gray: %2.2f (left target)', biasSizeX, biasSizeY,stabilizerGray);
        end
    end
    
    if (~isempty(rightGammaOut))
        % individual gamma point ratios (right target)
        for stabilizerGrayIndex = 1:stabilizerGrayLevelNum
            for biasSizeIndex = 1: biasSizesNum
                plot(squeeze(totalEnergy(stabilizerGrayIndex, biasSizeIndex,indices)), squeeze(ratiosRight(stabilizerGrayIndex, biasSizeIndex, indices)), 'ks-', 'LineWidth', 1, 'MarkerSize', 8, 'MarkerFaceColor', [1 1 1]); 
            end
        end
    end
    
    if (1==2)
    % The ratio of entire gamma curves
    for stabilizerGrayIndex = 1:stabilizerGrayLevelNum
        for biasSizeIndex = 1: biasSizesNum
            
            gammaCurveLeft = squeeze(leftGammaOut(stabilizerGrayIndex, biasSizeIndex,:));
            if (~isempty(rightGammaOut))
                gammaCurveRight = squeeze(rightGammaOut(stabilizerGrayIndex, biasSizeIndex,:)); 
            end
            
            % Entire curve ratio (minus first few points)
            curveRatiosLeft(stabilizerGrayIndex, biasSizeIndex) = 1.0/(gammaCurveLeft(indices) \ referenceGammaCurveLeft(indices));  
            
            condIndex = (stabilizerGrayIndex-1)* biasSizesNum + biasSizeIndex;
            lineColor = lineColors(condIndex,:);  
            plot(squeeze(totalEnergy(stabilizerGrayIndex, biasSizeIndex,indices)), repmat(squeeze(curveRatiosLeft(stabilizerGrayIndex, biasSizeIndex)), [1 size(totalEnergy,3)]), '-', 'LineWidth', 1, 'Color', lineColor);  
            
            if (~isempty(rightGammaOut))
                curveRatiosRight(stabilizerGrayIndex, biasSizeIndex) = 1.0/(gammaCurveRight(indices) \ referenceGammaCurveRight(indices));
                plot(squeeze(totalEnergy(stabilizerGrayIndex, biasSizeIndex,:)), repmat(squeeze(curveRatiosRight(stabilizerGrayIndex, biasSizeIndex)), [1 size(totalEnergy,3)]), 'k-', 'LineWidth', 1);
            end
       end
    end
    end
    
    
    hold off;
     
    xlabel('Panel activation (mean setting across all pixels)', 'FontName', 'Helvetica', 'FontSize', 10, 'FontWeight', 'bold');
    ylabel('luminance (cd/m2)', 'FontName', 'Helvetica', 'FontSize', 10, 'FontWeight', 'bold');
    ylabel('Gamma Ratio (measured/reference)', 'FontName', 'Helvetica', 'FontSize', 10, 'FontWeight', 'bold');
    set(gca, 'XLim', energyLims, 'YLim', [0.35 1.05], 'YTick', [0.1:0.1:1.1]);
    set(gca, 'FontName', 'Helvetica', 'FontSize', 8, 'Color', [1 1 1]);
    grid on;
    box on
    legend_handle = legend(legendMatrix{:}, 'Location', 'SouthWest');
    set(legend_handle, 'Box', 'on', 'FontName', 'Helvetica', 'FontSize', 12);
            
    %Print figure
    set(h2, 'Color', 0.9*[1 1 1]);
    set(h2,'PaperOrientation','Landscape');
    set(h2,'PaperUnits','normalized');
    set(h2,'PaperType', 'uslegal');
    set(h2, 'InvertHardCopy', 'off');
    set(h2,'PaperPosition', [0 0 1 1]);
    print(gcf, '-dpdf', 'Fig2.pdf');
    

    stabilizerGrayIndex = 2;
    biasSizeIndex = 3;
            
    condIndex = (stabilizerGrayIndex-1)* biasSizesNum + biasSizeIndex;
    lineColor = lineColors(condIndex,:);
            
    gammaCurveLeft = squeeze(leftGammaOut(stabilizerGrayIndex, biasSizeIndex,:));
            
    % Individual points ratio
    ratiosLeft = gammaCurveLeft./referenceGammaCurveLeft;
    % individual gamma point ratios (left target)
    
    h3 = figure(3);
    clf;
    
    subplot('Position', [0.07 0.68 0.91 0.3]);
    energy = squeeze(totalEnergy(stabilizerGrayIndex, biasSizeIndex,indices));
    dE = max(energy(:))*0.01;
    energyLims = [0.493 0.498];
    energyTicks = energyLims(1):0.001:energyLims(2);
    
    plot(squeeze(totalEnergy(stabilizerGrayIndex, biasSizeIndex,indices)),  ratiosLeft(indices), 'o-', 'MarkerSize', 8, 'MarkerFaceColor', lineColor, 'Color', lineColor*0.5, 'LineWidth', 1); 
    xlabel('');
    ylabel('Gamma Ratio (measured/reference)', 'FontName', 'Helvetica', 'FontSize', 10, 'FontWeight', 'bold');
    set(gca, 'XLim', energyLims, 'YLim', [0.35 1.05], 'YTick', [0.0:0.1:1.1], 'XTick', energyTicks);
    set(gca, 'FontName', 'Helvetica', 'FontSize', 8, 'Color', [1 1 1]);
    legend_handle = legend({'ratio'}, 'Location', 'SouthEast');
    set(legend_handle, 'Box', 'on', 'FontName', 'Helvetica', 'FontSize', 12);
    grid on;
    box on
    
    
    subplot('Position', [0.07 0.356 0.91 0.3]);
    hold on
    plot(squeeze(totalEnergy(stabilizerGrayIndex, biasSizeIndex,indices)),  gammaCurveLeft(indices)/max(gammaCurveLeft(indices)), 'o-', 'MarkerSize', 8, 'MarkerFaceColor', lineColor, 'Color', lineColor*0.5, 'LineWidth', 1); 
    plot(squeeze(totalEnergy(stabilizerGrayIndex, biasSizeIndex,indices)),  referenceGammaCurveLeft(indices)/max(referenceGammaCurveLeft(indices)), 'ks-', 'MarkerSize', 8,  'LineWidth', 1); 
    hold off
    xlabel('');
    ylabel('Normalized luminance', 'FontName', 'Helvetica', 'FontSize', 10, 'FontWeight', 'bold');
    set(gca, 'XLim', energyLims, 'YLim', [-0.05 1.05], 'YTick', [0.0:0.1:1.1], 'XTick', energyTicks);
    set(gca, 'FontName', 'Helvetica', 'FontSize', 8, 'Color', [1 1 1]);
    
    legend_handle = legend({'Measured (Bias region:420x420, Stabilizer gray: 0.33)', 'Reference(Bias region: none, Stabilizer gray: 0.0)'}, 'Location', 'SouthEast');
    set(legend_handle, 'Box', 'on', 'FontName', 'Helvetica', 'FontSize', 12);
    grid on;
    box on
    
    
    
    subplot('Position', [0.07 0.03 0.91 0.3]);
    hold on
    plot(squeeze(totalEnergy(stabilizerGrayIndex, biasSizeIndex,indices)),  gammaCurveLeft(indices)/max(gammaCurveLeft(indices)), 'o-', 'MarkerSize', 8, 'MarkerFaceColor', lineColor, 'Color', lineColor*0.5, 'LineWidth', 1); 
    plot(squeeze(totalEnergy(stabilizerGrayIndex, biasSizeIndex,indices)),  referenceGammaCurveLeft(indices)/max(referenceGammaCurveLeft(indices)), 'ks-', 'MarkerSize', 8,  'LineWidth', 1); 
    hold off
    set(gca, 'YScale', 'log', 'XScale', 'log');
    xlabel('Panel activation (mean setting across all pixels)', 'FontName', 'Helvetica', 'FontSize', 10, 'FontWeight', 'bold');
    ylabel('Normalized luminance', 'FontName', 'Helvetica', 'FontSize', 10, 'FontWeight', 'bold');
    set(gca, 'XLim', energyLims, 'YLim', [-0.05 1.05], 'YTick', [0.003 0.01 0.03 0.1 0.3 1.0], 'XTick', energyTicks);
    set(gca, 'FontName', 'Helvetica', 'FontSize', 8, 'Color', [1 1 1]);
    
    legend_handle = legend({'Measured (Bias region:420x420, Stabilizer gray: 0.33)', 'Reference(Bias region: none, Stabilizer gray: 0.0)'}, 'Location', 'SouthEast');
    
    set(legend_handle, 'Box', 'on', 'FontName', 'Helvetica', 'FontSize', 12);
    grid on;
    box on
    
    set(h3, 'Color', 0.9*[1 1 1]);
    set(h3,'PaperOrientation','Portrait');
    set(h3,'PaperUnits','normalized');
    set(h3,'PaperType', 'uslegal');
    set(h3, 'InvertHardCopy', 'off');
    set(h3,'PaperPosition', [0 0 1 1]);
    print(gcf, '-dpdf', 'Fig3.pdf');
    
    
    
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
    
    
    
    
    h1 = figure(99);
    clf;
    hold on;
    interpolatedSettingsValues = 0:0.001:1;
    c50 = [];
    if (exist('GammaTables.mat') == 2)
        s = load('GammaTables.mat');
        if isfield(s, 'data')
            data.gammaTables = s.data.gammaTables;
            data.c50 = s.data.c50;
        end
    else
       em = []; 
       save('GammaTables.mat', 'em')
    end
    for stabilizerGrayIndex = 1:stabilizerGrayLevelNum
        for biasSizeIndex = 1: biasSizesNum
            gammaCurveLeft = squeeze(gammaOutputLeft(stabilizerGrayIndex, biasSizeIndex, :));
            gammaCurveLeft = gammaCurveLeft / max(gammaCurveLeft(:));
            interpolatedGammaCurveLeft = interp1(gammaInputLeft, gammaCurveLeft, interpolatedSettingsValues, 'pchip');
            plot(gammaInputLeft, gammaCurveLeft, 'ks', 'MarkerSize', 12);
            plot(interpolatedSettingsValues, interpolatedGammaCurveLeft, 'r-');
            % search to find settings that produce 0.5 output
            [c,index] = min(abs(interpolatedGammaCurveLeft-0.5));
            plot(interpolatedSettingsValues(index), interpolatedGammaCurveLeft(index), 'bo');
            data.gammaTables(dataSetIndex,stabilizerGrayIndex,biasSizeIndex,:)= interpolatedGammaCurveLeft;
            data.c50(dataSetIndex,stabilizerGrayIndex,biasSizeIndex)       = interpolatedSettingsValues(index);
        end
    end
    drawnow;
    
    save('GammaTables.mat', 'interpolatedSettingsValues', 'data', '-append');
    whos('-file','GammaTables.mat')
    fprintf('saved data');
    pause;
    
    h1 = figure(1);
    figXo = 2560;   figYo = 360;
    figWidth = 700; figHeight = 860;
    set(h1, 'Position', [figXo figYo figWidth figHeight]);
    clf;
    
    lineColors = jet(stabilizerGrayLevelNum*biasSizesNum);
    
    % Subplot panel sizes and margins
    width = 0.88/(biasSizesNum);
    height = 0.74/(stabilizerGrayLevelNum);
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
            
            plot(gammaInputLeft, gammaCurveLeft, 'ko-', 'LineWidth', 1.0, 'MarkerSize', 8, 'MarkerFaceColor', lineColor, 'Color', lineColor*0.5);
            
            if (~isempty(gammaInputRight))
                gammaCurveRight = squeeze(gammaOutputRight(stabilizerGrayIndex, biasSizeIndex, :));
                if (maxGammaInputRight == minGammaInputRight)
                    gammaInput = gammaInputLeft;
                else
                    gammaInput = gammaInputRight;
                end
                
                hold on;
                plot(gammaInput, gammaCurveRight, 'ks-', 'LineWidth', 1.0, 'MarkerSize', 8, 'MarkerFaceColor', [1 1 1]);
                if (maxGammaInputRight == minGammaInputRight)
                    plot([minGammaInputRight, maxGammaInputRight], [0 mean(gammaCurveRight)], 'k-', 'LineWidth', 1.0);
                    plot(minGammaInputRight, 0, 'kv', 'LineWidth', 1.0, 'MarkerFaceColor', [0.6 0.6 0.6]);
                end
                hold off;
                legend_handle = legend({'Left', 'Right'}, 'FontName', 'Helvetica', 'FontSize', 6, 'Location', 'NorthWest');
                %set(legend_handle, 'Box', 'off')
            end
            
            set(gca, 'FontName', 'Helvetica', 'FontSize', 8, 'Color', [1 1 1]);
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
            
            if (stabilizerGrayIndex == stabilizerGrayLevelNum)
                xlabel('settings value', 'FontName', 'Helvetica', 'FontSize', 10, 'FontWeight', 'bold');
            else
                xlabel('');
            end
            
            
            title(sprintf('Stabilizer gray = %2.2f; \nBias WxH = %2.0fx%2.0f pxls.', stabilizerGray, biasSizeX, biasSizeY), 'FontName', 'Helvetica', 'FontSize', 8, 'BackgroundColor',[.99 .99 .48], 'EdgeColor', [0 0 0]);
           
        end
        
        ylabel('');
        set(gca, 'YTickLabel', []);
        set(gca, 'YLim', [0 maxGammaOutputLeft]);
        set(gca, 'XTick', [0:0.2:1.0], 'YTick', [0:100:1000]);
        set(gca, 'FontName', 'Helvetica', 'FontSize', 8, 'Color', [1 1 1]);
        grid on;
        box on

    end
    
   
        
    
    % Print figure
    set(h1, 'Color', 0.9*[1 1 1]);
    set(h1,'PaperOrientation','Portrait');
    set(h1,'PaperUnits','normalized');
    set(h1, 'InvertHardCopy', 'off');
    set(h1,'PaperPosition', [0 0 1 1]);
    print(gcf, '-dpdf', 'Fig1.pdf');
    
    
end


