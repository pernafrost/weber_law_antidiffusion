%% Movement response of a particle belonging to a Gaussian concentration
% This script generates the x positions of a particle that belongs to a
% stable Gaussian shaped group. The x positions are autocorrelated in time
% The number of simulated trajectories is nTrajectories and the length of
% each trajectory is tL.



clear all
close all

rng(753)

%% parameters for simulation
corLen = 10; % autocorrelation length
tL = 1024; % trajectory length
lambda = 1; % this is analogous to the sigma of the gaussian distribution
mu = 0; % this is the mean of the gaussian distribution
nTrajectories = 50; % number of simulated trajectories



%% parameters for visualization
useBlackBackground = 0;
figureFolder = 'Figures';
isSavingFigures = 1;
markerSize = 10;
markerSizeSmall = 8;
useFontSize = 20;
prefixToFileNames = '';

if useBlackBackground == 1
    colordef('black');
    figureBackgroundColour = [0, 0, 0];
    markerFaceColour = [1, 0.5, 0.3];
    lineColour = [1, 1, 1];
    currentColourMap = hot(1000);
    currentColourMap = currentColourMap(end:-1:1, :);
    lineWidth = 2;
else
    colordef('white');
    figureBackgroundColour = [1, 1, 1];
    markerFaceColour = [0.5, 0.5, 0.5];
    lineColour = [0, 0, 0];
    currentColourMap = hot(1000);
    lineWidth = 2;
end


if ~exist(figureFolder, 'dir') && isSavingFigures == 1
    mkdir(figureFolder);
end


%% correlation matrix
[ii,jj] = meshgrid(1:tL);
K = lambda^2*exp(-(ii-jj).^2/corLen^2);
% the condition below is only necessary if we want the end point to be the
% same as the first point
C = K; %- K([1, end], :)* inv([K(1,[1, end]); K(end,[1 ,end])]) * K(:, [1, end]);


%% The theoretical gaussian distribution that describes the group density
xBins = -3:0.2:3;
centerXBins = xBins(1:end-1) + diff(xBins)/2;
xBinsSmall = -3:0.01:3;
centerXBinsSmall = xBinsSmall(1:end-1) + diff(xBinsSmall)/2;
norm1 = normpdf(centerXBinsSmall,mu,lambda);
norm1Grad = diff(normpdf(xBinsSmall,mu,lambda)) ./ diff(xBinsSmall);

% I expect that for xBins small enough the two columns are very similar
% [centerXBinsSmall'/lambda^2, -(norm1Grad ./ norm1)']


allXPositions = nan(nTrajectories * (tL-1), 1); % I exclude the last position because I don't know the displacement corresponding to that position
allDeltaX = nan(nTrajectories * (tL-1), 1);
%% generate random trajectories
for tt = 1:nTrajectories
    rn = mvnrnd(mu * ones(1,tL), C);
    % the movement of the particle at each time step is
    deltaX = diff(rn);
    
    allXPositions((tt - 1) * (tL - 1) + 1: tt * (tL - 1)) = rn(1:end-1);
    allDeltaX((tt - 1) * (tL - 1) + 1: tt * (tL - 1)) = deltaX;
    
    % this is a check: muHat and sigmaHat should be equal to mu and sigma
    % [muHat,sigmaHat] = normfit(rn)
end

[nInBin, valXInBin] = histc(allXPositions, xBins);

nXBins = length(xBins) - 1;
nCountsX = NaN(nXBins, 1);
meanDeltaX = NaN(nXBins, 1);
stdDeltaX = NaN(nXBins, 1);
for bb=1:nXBins
    valX = find(valXInBin == bb);
    nCountsX(bb) = length(valX);
    meanDeltaX(bb) = mean(allDeltaX(valX));
    stdDeltaX(bb) = std(allDeltaX(valX)) / sqrt(nCountsX(bb) - 1);
end


regressionSlope = allXPositions \ allDeltaX; % allDeltaX \ allXPositions;
% regressionSlope2 = polyfit(allXPositions, allDeltaX, 1)

% If there is no correlation (corLen = 1) then at the next time step we
% should have deltaX = -x on average
% empirically deltaX = - x * 2 / pi

%% Plot (average) delta X vs position
myFig = figure;
set(myFig, 'Position', [1 1 1000 600]);
set(myFig, 'Color', figureBackgroundColour);
hold on;
plot(centerXBins, regressionSlope * centerXBins, '-', 'Color', markerFaceColour, 'LineWidth', lineWidth);

errorbar(centerXBins, meanDeltaX, stdDeltaX, 'o', 'Color', markerFaceColour);
plot(centerXBins, meanDeltaX, 'Marker', 'o', 'MarkerFaceColor', markerFaceColour, 'Color', markerFaceColour, 'LineStyle', 'none', 'MarkerSize', markerSize, 'MarkerEdgeColor', lineColour, 'LineWidth', lineWidth);


set(gca, 'LineWidth', lineWidth, 'FontSize', useFontSize, 'TickDir', 'out', 'FontName', 'Arial');
set(gca, 'XLim', [xBins(1), xBins(end)]);
set(gca, 'YLim', [-0.06, 0.06]);
xlabel('x'); ylabel('\Delta x');
text(2.2, 0.05, '(c)', 'FontName', 'Arial', 'FontSize', 50);
axis square;

if isSavingFigures
    set(gcf,'PaperPositionMode','auto');
    figureFileName = fullfile(figureFolder, [prefixToFileNames, 'delta_x_vs_x.eps']);
    print(gcf, '-depsc2', '-tiff', '-r300', figureFileName);
end




%% plot one trajectory
myFig2 = figure;
set(myFig2, 'Position', [1 1 1200 560]);
set(myFig2, 'Color', figureBackgroundColour);
subplot(1,8,1:3);
hold on;
plot(rn(1:tL/2), 'Marker', 'o', 'MarkerFaceColor', markerFaceColour, 'Color', markerFaceColour, 'LineStyle', '-', 'MarkerSize', markerSizeSmall, 'MarkerEdgeColor', lineColour, 'LineWidth', lineWidth);
set(gca, 'LineWidth', lineWidth, 'FontSize', useFontSize, 'TickDir', 'out', 'FontName', 'Arial');
set(gca, 'YLim', [xBins(1), xBins(end)]);
set(gca, 'XLim', [-20, tL/2 + 20 + 1]);
xlabel('time t'); ylabel('x position');
text(10, 2.5, '(a)', 'FontName', 'Arial', 'FontSize', 50);
subplot(1,8,4.3)
set(gca, 'LineWidth', lineWidth, 'FontSize', useFontSize, 'TickDir', 'out', 'FontName', 'Arial');
plot(centerXBinsSmall, norm1, '-', 'Color', lineColour, 'LineWidth', lineWidth);
set(gca, 'XLim', [xBins(1), xBins(end)]);
set(gca, 'LineWidth', lineWidth, 'FontSize', useFontSize, 'TickDir', 'out', 'FontName', 'Arial');
xlabel('x position'); ylabel('prob. density');
box off;
view(90,90) % swap the x and y axis
text(-2.5, 0.13, '(b)', 'FontName', 'Arial', 'FontSize', 50);


% subplot with regression
subplot(1,8,6:8)
plot(centerXBins, regressionSlope * centerXBins, '-', 'Color', markerFaceColour, 'LineWidth', lineWidth);
hold on;
errorbar(centerXBins, meanDeltaX, stdDeltaX, 'o', 'Color', markerFaceColour);
plot(centerXBins, meanDeltaX, 'Marker', 'o', 'MarkerFaceColor', markerFaceColour, 'Color', markerFaceColour, 'LineStyle', 'none', 'MarkerSize', markerSize, 'MarkerEdgeColor', lineColour, 'LineWidth', lineWidth);

box off;
set(gca, 'LineWidth', lineWidth, 'FontSize', useFontSize, 'TickDir', 'out', 'FontName', 'Arial');
set(gca, 'XLim', [xBins(1), xBins(end)]);
set(gca, 'YLim', [-0.06, 0.06]);
set(gca, 'XTick', [-2, -1, 0, 1, 2]);
xlabel('x position'); ylabel('\Delta x');
text(-1.8, 0.05, '(c)', 'FontName', 'Arial', 'FontSize', 50);
% axis square;



if isSavingFigures
    set(gcf,'PaperPositionMode','auto');
    figureFileName = fullfile(figureFolder, [prefixToFileNames, 'correlated_trajectory_in_gaussian_group.eps']);
    print(gcf, '-depsc2', '-tiff', '-r300', figureFileName);
end



warndlg('Now running a loop on different values of lambda and corLen', 'This will take a lot of time');

%% create heatmap with all values of slope vs. lambda and corLen

lambdaVals = 0.5:0.5:15;
corLenVals = 1:1:20;
allRegressionSlopes = nan(length(lambdaVals), length(corLenVals));
for ll = 1:length(lambdaVals)
    lambda = lambdaVals(ll);
    for cc = 1:length(corLenVals)
        corLen = corLenVals(cc);
        fprintf('iterating at lambda=%f and corLen=%f\n', lambda, corLen);
        
        % correlation matrix
        [ii,jj] = meshgrid(1:tL);
        K = lambda^2*exp(-(ii-jj).^2/corLen^2);
        % the condition below is only necessary if we want the end point to be the
        % same as the first point
        C = K; %- K([1, end], :)* inv([K(1,[1, end]); K(end,[1 ,end])]) * K(:, [1, end]);
        
        
        
        allXPositions = nan(nTrajectories * (tL-1), 1); % I exclude the last position because I don't know the displacement corresponding to that position
        allDeltaX = nan(nTrajectories * (tL-1), 1);
        for tt = 1:nTrajectories
            rn = mvnrnd(mu * ones(1,tL), C);
            % the movement of the particle at each time step is
            deltaX = diff(rn);
            
            allXPositions((tt - 1) * (tL - 1) + 1: tt * (tL - 1)) = rn(1:end-1);
            allDeltaX((tt - 1) * (tL - 1) + 1: tt * (tL - 1)) = deltaX;
            
            % this is a check: muHat and sigmaHat should be equal to mu and sigma
            % [muHat,sigmaHat] = normfit(rn)
        end
        regressionSlope(ll, cc) = allXPositions \ allDeltaX;
        
    end
end

% save('regressionSlope.mat', regressionSlope, lambdaVals, corLenVals);

% with imagesc the ticks are at the right position
figure, pcolor(corLenVals, lambdaVals, regressionSlope); colorbar();
shading flat;
ylabel('lambda'); xlabel('corr. len.');


