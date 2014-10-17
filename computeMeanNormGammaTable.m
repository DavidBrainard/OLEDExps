function computeMeanNormGammaTable()
    s = load('GammaTables.mat');
    s.data
    gammaTables = s.data.gammaTables;
    gammaTables2 = reshape(gammaTables, [size(gammaTables,1)*size(gammaTables,2)*size(gammaTables,3) size(gammaTables,4)]);
    meanGammaCurve = squeeze(mean(gammaTables2,1));
    
    c50 = s.data.c50;
    c50 = c50(:);
    
    h1 = figure(99);
    set(h1, 'Position', [300 300 870 435]);
    clf;
    subplot(1,2,1);
    hold on;
    for i = 1:size(gammaTables,1)
        for j = 1:size(gammaTables,2)
            for k = 1:size(gammaTables,3)
                plot(s.interpolatedSettingsValues, squeeze(gammaTables(i,j,k,:)), 'k-');
            end
        end
    end
    plot(s.interpolatedSettingsValues, meanGammaCurve, 'r-', 'LineWidth', 2);
    drawnow;
    set(gca, 'FontName', 'Helvetica', 'FontSize', 14, 'Color', [1 1 1]);
    set(gca, 'XTick', [0:0.1:1], 'YTick', [0:0.1:1]);
    xlabel('settings value', 'FontName', 'Helvetica', 'FontSize', 16, 'FontWeight', 'bold');
    ylabel('normalized output', 'FontName', 'Helvetica', 'FontSize', 16, 'FontWeight', 'bold');
    title(sprintf('Gamma curve population (%d) & mean', size(gammaTables2,1)), 'FontName', 'Helvetica', 'FontSize', 14, 'FontWeight', 'bold');
    grid on;
    box on
    axis 'square'
    
    
    subplot(1,2,2);
    binPositions = 0.71:0.002:0.77;
    h = histogram(c50,binPositions);
    h.FaceColor = [0 0.5 0.5];
    h.EdgeColor = 'r';
    
    set(gca, 'FontName', 'Helvetica', 'FontSize', 14, 'Color', [1 1 1]);
    set(gca, 'XTick', [0.7:0.01:0.76], 'XLim', [0.72 0.75]);
    xlabel('c50', 'FontName', 'Helvetica', 'FontSize', 16, 'FontWeight', 'bold');
    ylabel('Count', 'FontName', 'Helvetica', 'FontSize', 16, 'FontWeight', 'bold');
    title(sprintf('Distribution of c50 (median:%2.3f)', median(c50)), 'FontName', 'Helvetica', 'FontSize', 14, 'FontWeight', 'bold');
    grid on;
    box on
    axis 'square'
    
    % Print figure
    set(h1, 'Color', 0.9*[1 1 1]);
    set(h1,'PaperOrientation','Landscape');
    set(h1,'PaperUnits','normalized');
    set(h1, 'InvertHardCopy', 'off');
    set(h1,'PaperPosition', [0 0 1 1]);
    print(gcf, '-dpdf', 'PreliminaryGammaCurvePopulation.pdf');
    
    gammaFunction.input = s.interpolatedSettingsValues;
    gammaFunction.output = meanGammaCurve;
    save('GammaFunction.mat', 'gammaFunction');
    
end

