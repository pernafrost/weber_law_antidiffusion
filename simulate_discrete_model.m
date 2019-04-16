

clear all;
close all

useBlackBackground = 0;
isSavingFigures = 0;
figureFolder = 'figures';
len = 512;
% im = 1000 * ones(len,1);

if ~exist(figureFolder, 'dir') && isSavingFigures == 1
    mkdir(figureFolder);
end


% im(15) = 701;
nTestedCycles = [3, 6, 12, 24];% 1:1:50;
testedCycle = nan(size(nTestedCycles));
localMax = nan(size(nTestedCycles));
for jj = 1:length(nTestedCycles)
    nCycles = nTestedCycles(jj);
    testedCycle(jj) = nCycles;
    
    im= 100+1*sin((1:len)*nCycles*2*pi/len);
    
    D = 0.01;
    
    gamma = 0.02;%D/ (1- 3*D);
    
    L = [1, -2, 1];
    
    if useBlackBackground == 1
        colordef('black');
    end
    
    for tt = 1:100
        im = wextend('1D','ppd',im, 1);
        im = im + conv(im, L, 'same')*D;
        im = im(2:end-1);
        
        im = wextend('1D','ppd',im, 1);
        im = im - conv(im, L, 'same')*gamma;
        im = im(2:end-1);
        
        sum(im(:))
        
        plot(im, 'LineWidth', 2);
        % set(gca, 'YLim', [0, 2000])
        box off;
        set(gca, 'TickDir', 'out', 'LineWidth', 2);
        set(gca, 'FontName', 'Arial', 'FontSize', 18);
        set(gcf, 'InvertHardCopy', 'off');
        set(gca, 'XLim', [1, len])
        drawnow;
        
        if isSavingFigures
            figureFileName = fullfile(figureFolder, sprintf('sine_wave%05d_cycles_t%d.png', nCycles, tt));
            print(gcf, '-dpng', '-r300', figureFileName);
        end
        
    end
    
    localMax(jj) = max(im);
end



figure, loglog(1./testedCycle, localMax/1000, '.-')


