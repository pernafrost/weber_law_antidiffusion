

clearvars;
close all;


% nVisibleParticles = 11; % these are only for illustrative purpose
D = 0.3;
useBlackBackground = 0;
figureFolder = 'Figures';
isSavingFigures = 1;
markerSize = 16;
lineWidth = 4;
markerLineWidth = 2;
useFontSize = 36;

nParticles = 24;
nBins = 2;
for isWeber = 0 % if use Weber Law type response instead of diffusion
    % isWeber == 0: only diffusion; isWeber > 0 : Weber Law up gradient;
    % isWeber < 0: Weber Law down gradient
    
    if isSavingFigures && ~exist(figureFolder, 'dir')
        mkdir(figureFolder);
    end
    
    if useBlackBackground == 1
        colordef('black');
        figureBackgroundColour = [0, 0, 0];
        markerFaceColour = [1, 0.5, 0.3];
        markerFacePossibleColours = jet(nBins);
        % rng(24);
        markerFacePossibleColours = markerFacePossibleColours(randperm(nBins), :);
        lineColour = [1, 1, 1];
        currentColourMap = hot(1000);
        currentColourMap = currentColourMap(end:-1:1, :);
        % visibleParticleColours = ones(nVisibleParticles, 3); % jet(nVisibleParticles);
        % visibleParticleColours(1,:) = [1, 0.3, 0.3];
        markerFacePossibleColours(1,:) = [0.3, 0.3, 0.3];
        markerFacePossibleColours(4,:) = [0.3, 0.3, 0.3];
        markerFacePossibleColours(5,:) = [0.3, 0.3, 0.3];
        
    else
        colordef('white');
        figureBackgroundColour = [1, 1, 1];
        markerFaceColour = [0.85, 0.85, 0.85];
        markerFacePossibleColours = repmat(markerFaceColour,nBins,1);
        lineColour = [0, 0, 0];
        currentColourMap = hot(1000);
        % visibleParticleColours = zeros(nVisibleParticles, 3); % jet(nVisibleParticles); %
        % visibleParticleColours(1,:) = [0.8, 0, 0];
        markerFacePossibleColours(2,:) = [0.8, 0.2, 0];
        markerFacePossibleColours(3,:) = [0, 0.1, 0.6];
    end
    
    
    
    % rng(2); % rng(1)
    xBins = 0:1:(nBins + 1);
    xCenterBins = xBins(1:end-1) + diff(xBins)/2;
    xTickLabels = {'l', 'r'};
    
    xLims = [0, nBins]; xLimsExtended = xLims + [-0.02, +0.02];
    yLims = [0, 1]; yLimsExtended = yLims + [-0.02, 0];
    yBarHeight = 0.6;
    
    %% make a drawing of the discrete diffusion process
    myFig = figure;
    set(myFig, 'Position', [1 1 1000 600]);
    set(myFig, 'Color', figureBackgroundColour);
    hold on;
    set(gca, 'TickDir', 'out', 'LineWidth', 1.5, 'box', 'off');
    set(gca, 'FontName', 'arial', 'FontSize', useFontSize);
    hxl = xlabel('Spatial position x');
    set(gca, 'Clipping', 'off');
    
    
    % offset xlabel
    offset = 0.05;
    set(hxl, 'Units', 'Normalized');
    pos = get(hxl, 'Position');
    set(hxl, 'Position', pos + [0, -offset, 0]);
    
    % offset axes
    offset = 0.15;
    pos = get(gca, 'Position');
    set(gca, 'Position', pos + [0, offset, 0, 0]);
    
    
    % ylabel('');
    stem(xBins, yBarHeight * ones(size(xBins)), 'Marker', 'none', 'LineWidth', lineWidth, 'Color', lineColour);
    plot(xLims, [0,0], 'k-', 'LineWidth', lineWidth)
    set(gca, 'YTick', []);
    xticks(xCenterBins(1:end-1),xTickLabels,'rotation',0,'fontsize',useFontSize)
    % set(gca, 'XTick', xCenterBins, 'XTickLabel', xTickLabels);
    set(gca, 'XLim', xLimsExtended, 'YLim', yLimsExtended);
    
    set(gca, 'visible', 'off');
    ah = gca;
    set(findall(gca, 'type', 'text'), 'visible', 'on')
    % set(findall(gca, 'type', 'Tick'), 'visible', 'on')
    % gca.XTick = 'visible';
    
    % plot particle positions
    
    particleRadius = 0.13; % I don't know the particle radius corresponding to MarkerSize, but I make up another radius for plotting the particles
    if particleRadius > 0.5 % I cannot have particles as large as the entire bin
        particleRadius = 0.45;
    end
    
    % xBinRand = randi(nBins, [nParticles, 1]);
    
    % put more particles in the central bin
    % xBinRand(1:14) = 3; xBinRand(4) = 2;
    
    xBinRand = [ones(6,1); 2*ones(18,1)];
    
    markerFaceColours = markerFacePossibleColours(xBinRand, :);
    
    
    xPosRandInBin = rand(nParticles,1) * (1-2*particleRadius) + particleRadius;
    xPos = xBinRand - 1 + xPosRandInBin;
    yPos = rand(nParticles, 1) * yBarHeight * 5/6;
    
    floorRadius = 0.04;
    
    figureCounter = 1;
    % update positions
    borderAvoidanceYesNo = 1;
    for kk = 1:nParticles*12005
        
        updateParticleID = randi(nParticles);
        
        % calculate forces
        dx = (xPos - xPos(updateParticleID));
        dy = (yPos - yPos(updateParticleID));
        d = (dx.^2 + dy.^2).^0.5; % squared distance
        
        avoidancex = 0;
        avoidancey = 0;
        interacting = setdiff(find(d <= particleRadius), updateParticleID); % interacting particles are all within particleRadius, except the particle itself
        closestNeighbour = NaN;
        if ~isempty(interacting)
            [~, c] = min(d(interacting));
            closestNeighbour = interacting(c);
            avoidancex = -dx(closestNeighbour);
            avoidancey = -dy(closestNeighbour);
        end
        distFromClosestBorder = xPos(updateParticleID) - round(xPos(updateParticleID));
        borderAvoidancex = 0;
        if abs(distFromClosestBorder) <= particleRadius && yPos(updateParticleID) <= yBarHeight + particleRadius ...
                || xPos(updateParticleID) <= particleRadius || xPos(updateParticleID) >= 5 - particleRadius
            borderAvoidancex = sign(distFromClosestBorder);
        end
        
        floorAvoidancey = 0;
        if yPos(updateParticleID) <= floorRadius
            floorAvoidancey = 1;
        end
        
        if yPos(updateParticleID) <= yBarHeight + particleRadius
            noisex = randn(1);
        else
            noisex = 2 * randn(1);
        end
        noisey = 2 * randn(1);
        gravity = -rand(1);
        
        correctiony = 0;
        if yPos(updateParticleID) > yBarHeight * 5/6
            correctiony = -1.5;
        end
        if yPos(updateParticleID) < floorAvoidancey
            correctiony = 1.5;
        end
        correctionx = 0;
        if xPos(updateParticleID) > nBins - particleRadius
            correctionx = -1.5;
        end
        if xPos(updateParticleID) < 0 + particleRadius
            correctionx = +1.5;
        end
        
        particleMotionx = 0.01 * borderAvoidanceYesNo * avoidancex + 0.01 * borderAvoidanceYesNo * borderAvoidancex + 0.005 * noisex + 0.01 * correctionx;
        particleMotiony = 0.01 * borderAvoidanceYesNo * avoidancey + 0.0001 * gravity + 0.01*floorAvoidancey + 0.005 * noisey + 0.01 * correctiony;
        
        % particleMotionx = 0.00001 * borderAvoidancex + 0.002 * noisex + 0.02 * correctionx;
        % particleMotiony = 0.01 * avoidancey + 0.0002 * gravity + 0.01*floorAvoidancey + 0.002 * noisey + 0.02 * correctiony;
        
        
        xPos(updateParticleID) = xPos(updateParticleID) + particleMotionx;
        yPos(updateParticleID) = yPos(updateParticleID) + particleMotiony;
        
        %     if xPos(updateParticleID) < 0
        %         xPos(updateParticleID) = -xPos(updateParticleID);
        %     end
        %     if yPos(updateParticleID) < 0
        %         yPos(updateParticleID) = - yPos(updateParticleID);
        %     end
        %
        %     if xPos(updateParticleID) > nBins
        %         xPos(updateParticleID) = xPos(updateParticleID) - nBins;
        %     end
        %     if yPos(updateParticleID) > 1
        %         yPos(updateParticleID) = yPos(updateParticleID) - 1;
        %     end
        
        nInEachBin = histc(xPos, 0:1:nBins);
        nTot = sum(nInEachBin);
        
        if mod(kk, nParticles*6000) == 0 % || mod(kk, nParticles*6000) == nParticles*1000
            kk
            cla;
            plot(xLims, [0,0], '-', 'LineWidth', lineWidth, 'Color', lineColour);

            stem(xBins, 0.65 * ones(size(xBins)), 'Marker', 'none', 'LineWidth', lineWidth, 'Color', lineColour);
            %        plot(xPos, yPos, 'Marker', 'o', 'Color', lineColour, 'MarkerFaceColor', markerFaceColour, 'MarkerSize', markerSize, 'LineWidth', markerLineWidth, 'LineStyle', 'none');
            
            for jj = 1: nParticles
                plot(xPos(jj), yPos(jj), 'Marker', 'o', 'Color', lineColour, 'MarkerFaceColor', markerFaceColours(jj,:), 'MarkerSize', markerSize, 'LineWidth', markerLineWidth, 'LineStyle', 'none');
            end
            % plot(xPos(updateParticleID), yPos(updateParticleID), 'Marker', 'o', 'Color', 'r', 'MarkerFaceColor', 'r', 'MarkerSize', markerSize, 'LineWidth', markerLineWidth, 'LineStyle', 'none');
            % plot(xPos(interacting), yPos(interacting), 'Marker', 'o', 'Color', 'g', 'MarkerFaceColor', 'g', 'MarkerSize', markerSize, 'LineWidth', markerLineWidth, 'LineStyle', 'none');
            % if isfinite(closestNeighbour)
            %     plot(xPos(closestNeighbour), yPos(closestNeighbour), 'Marker', 'o', 'Color', 'k', 'MarkerFaceColor', 'g', 'MarkerSize', markerSize, 'LineWidth', markerLineWidth, 'LineStyle', 'none');
            % end
            % quiver(xPos(updateParticleID), yPos(updateParticleID), avoidancex, avoidancey);
            xticks(xCenterBins(1:end-1), xTickLabels, 'rotation', 0, 'fontsize', useFontSize)
            
%             for jj = 1: (length(nInEachBin) - 1)
%                 text(xCenterBins(jj) - 0.1, 0.7, sprintf('N = %d', nInEachBin(jj)), 'FontName', 'arial', 'FontSize', useFontSize);
%             end
            
            axis tight
            set(gca, 'XLim', xLimsExtended, 'YLim', yLimsExtended);
            drawnow;
            
            if isSavingFigures && borderAvoidanceYesNo
                figureFileName = sprintf('%s%cillustrative_drawing%05d.eps', figureFolder, filesep, figureCounter);
                if isWeber < 0
                    figureFileName = sprintf('%s%cillustrative_drawing_down_gradient%05d.eps', figureFolder, filesep, figureCounter);
                end
                if isWeber > 0
                    figureFileName = sprintf('%s%cillustrative_drawing_up_gradient%05d.eps', figureFolder, filesep, figureCounter);
                end
                
                set(gcf, 'InvertHardCopy', 'off');
                set(gcf, 'PaperPositionMode', 'auto');
                print(gcf, '-depsc2', '-tiff', '-r300', figureFileName); print(gcf, '-dpng', '-r300', [figureFileName(1:end-3), 'png'])
                figureCounter = figureCounter + 1;
                % rng(2)
            end
            
           
            
            if sum(nInEachBin) < nParticles
                break;
            end
            
            % borderAvoidanceYesNo = mod(borderAvoidanceYesNo + 1, 2)
            
            % pause(1)
        end
        
    end
    
    cla;
    stem(xBins, 1 * ones(size(xBins)), 'Marker', 'none', 'LineWidth', lineWidth, 'Color', lineColour);
    plot(xPos, yPos, 'Marker', 'o', 'Color', lineColour, 'MarkerFaceColor', markerFaceColour, 'MarkerSize', markerSize, 'LineWidth', markerLineWidth, 'LineStyle', 'none');
    xticks(xCenterBins(1:end-1), xTickLabels, 'rotation', 0, 'fontsize', useFontSize)
    
%     for jj = 1: (length(nInEachBin) - 1)
%         text(xCenterBins(jj)-0.2, 0.94, sprintf('N = %d', nInEachBin(jj)), 'FontName', 'arial', 'FontSize', useFontSize);
%     end
    
    drawnow;
    
    
    close all;
end
% if a particle is on top of the grid, assign a new position to it
% if a particle is too close to a different particle, then move it


%
% if isSavingFigures
%     figureFileName = sprintf('%s%cillustrative_drawing.eps', figureFolder, filesep);
%     set(gcf, 'InvertHardCopy', 'off');
%     set(gcf, 'PaperPositionMode', 'auto');
%     print(gcf, '-depsc2', '-tiff', '-r300', figureFileName); print(gcf, '-dpng', '-r300', [figureFileName(1:end-3), 'png'])
% end
%


