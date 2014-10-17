function plotC50sFromPreliminaryExperiments
    c50(1,:) = [0.7310    0.7360    0.7440    0.7340    0.7410    0.7460    0.7390    0.7430    0.7470    0.7410    0.7450    0.7470];
    c50(2,:) = [0.7310    0.7360    0.7410    0.7330    0.7380    0.7430    0.7370    0.7430    0.7440    0.7400    0.7450    0.7460];
    c50(3,:) = [0.7290    0.7350    0.7380    0.7320    0.7380    0.7420    0.7360    0.7400    0.7430    0.7400    0.7430    0.7430];
    c50(4,:) = [0.7250    0.7330    0.7390    0.7270    0.7350    0.7420    0.7340    0.7400    0.7440    0.7360    0.7430    0.7460];
    c50(5,:) = [0.7260    0.7320    0.7400    0.7250    0.7330    0.7400    0.7270    0.7350    0.7410    0.7290    0.7370    0.7420];
    c50(6,:) = [0.7340    0.7390    0.7440    0.7340    0.7410    0.7440    0.7370    0.7410    0.7430    0.7350    0.7410    0.7440];
    c50(7,:) = [0.7330    0.7380    0.7440    0.7340    0.7400    0.7450    0.7360    0.7400    0.7440    0.7340    0.7400    0.7450];
    c50(8,:) = [0.7340    0.7390    0.7430    0.7340    0.7390    0.7440    0.7360    0.7410    0.7410    0.7320    0.7390    0.7440];
    
    c50 = c50(:);
    h1 = figure(10);
    binPositions = 0.71:0.002:0.77;
    counts = hist(c50,binPositions);
    h = findobj(gca,'Type','patch');
    h.FaceColor = [0.95 0.2 0.5];
    h.EdgeColor = 'w';

    set(gca, 'FontName', 'Helvetica', 'FontSize', 14, 'Color', [1 1 1]);
    set(gca, 'XTick', [0.7:0.01:0.76]);
    xlabel('c50', 'FontName', 'Helvetica', 'FontSize', 16, 'FontWeight', 'bold');
    ylabel('Count', 'FontName', 'Helvetica', 'FontSize', 16, 'FontWeight', 'bold');
    title(sprintf('Distribution of c50 (median:%2.3f)', median(c50)), 'FontName', 'Helvetica', 'FontSize', 14, 'FontWeight', 'bold');
    grid on;
    box on
    
    % Print figure
    set(h1, 'Color', 0.9*[1 1 1]);
    set(h1,'PaperOrientation','Portrait');
    set(h1,'PaperUnits','normalized');
    set(h1, 'InvertHardCopy', 'off');
    set(h1,'PaperPosition', [0 0 1 1]);
    print(gcf, '-dpdf', 'c50.pdf');
    
    
end

