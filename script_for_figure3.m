% This script allows to reproduce the analysis shown in figure 3.
% The script reads an initial image (that can be considered as a
% two-dimensional density landscape) and then applies diffusion and
% anti-diffusion processes to it.
% There are three main possible configurations


clear all;
close all;

rng(10000); % I use the same seed
% rng shuffle; % I use a different seed
isSavingFigures = 1;

configuration = 'Dora'; % 'circles';% 'Dora'; 'cameraman'; ''
switch configuration
    
    case {'Dora', 'cameraman', 'Parkstead'}
        invertColours = 0; % if true invert the colours of the image before starting
        nVisibleParticles = 50; % these are particles superposed to the image to
        % visualize the diffusion and anti-diffusion process: the image does not
        % interact with these particles, but they do interact with the image.
        
        onlyPlotAtTheEndOfCycle = 0;
        
        colourMapForParticles = jet(nVisibleParticles);
        
        imTitle = [configuration, '_diff']; % title of saved figures
        imCounter = 1;
        figureFolder = 'figures';
        exportImageResolution = '-r300';
        isPlottingTrajectories = 1; clearTrajectoriesBeforeAntidiffusion = 1;
        
        iMax = 255; % the maximum value in the image at the beginning
        isBlurringWithLaplacian = 1; % using the Laplacian convolution is equivalent to taking a mean field approach
        isDeBlurringWithLaplacian = 0; % when deblurring, the laplacian convolution can produce images with negative values
        useMeanFieldForDeblurring = 1; % there isn't a useMeanFieldForBlurring because this would be equivalent to the Laplacian; in the case of deblurring
        % the further difference is whether the single particles that move
        % according to Weber's law do it probabilistically, or on average. If
        % isDeBlurringWithLaplacian, this useMeanField is ignored
%         reactionTime = 0; % if the reaction time is greater than zero, the
%         % de-blurring particles sense their environment at time t -
%         % reactionTime and move according to what they have sensed at time
%         % t
        avoidSaturationWhenDeblurring = 1; % if avoid saturation, then the particles won't move to a cell with more than iMax particles
        nBlurringSteps = 300; % 300
        nDeBlurringSteps = 300;
        nBlurringDeblurringCycles = 1;
        
        D = 0.01;% 0.01; % diffusion coefficient
        gamma = D ./ (1 -5*D); % anti-diffusion coefficient
        % gamma = D;% * sqrt(2);
        
        
        try
            imFileName = [configuration, '.png']; %uigetfile(); % 'cameraman.tif';%'../figures/File0493.png';%
        catch
            imFileName = [configuration, '.tif'];
        end
        % load image
        im = imread(imFileName);
        
        
        
    case 'circles'
        invertColours = 0; % if true invert the colours of the image before starting
        nVisibleParticles = 50; % these are particles superposed to the image to
        % visualize the diffusion and anti-diffusion process: the image does not
        % interact with these particles, but they do interact with the image.
        
        onlyPlotAtTheEndOfCycle = 0;
        
        
        colourMapForParticles = jet(nVisibleParticles);
        
        imTitle = 'circles_diff'; % title of saved figures
        imCounter = 1;
        figureFolder = 'figures';
        exportImageResolution = '-r300';
        isPlottingTrajectories = 1; clearTrajectoriesBeforeAntidiffusion = 1;
        
        iMax = 255; % the maximum value in the image at the beginning
        isBlurringWithLaplacian = 1; % using the Laplacian convolution is equivalent to taking a mean field approach
        isDeBlurringWithLaplacian = 0; % when deblurring, the laplacian convolution can produce images with negative values
        useMeanFieldForDeblurring = 1; % there isn't a useMeanFieldForBlurring because this would be equivalent to the Laplacian; in the case of deblurring
        % the further difference is whether the single particles that move
        % according to Weber's law do it probabilistically, or on average. If
        % isDeBlurringWithLaplacian, this useMeanField is ignored
        
%         reactionTime = 0; % if the reaction time is greater than zero, the
%         % de-blurring particles sense their environment at time t -
%         % reactionTime and move according to what they have sensed at time
%         % t
        avoidSaturationWhenDeblurring = 1; % if avoid saturation, then the particles won't move to a cell with more than iMax particles
        nBlurringSteps = 300; % 300
        nDeBlurringSteps = 300;
        nBlurringDeblurringCycles = 1;
        
        D = 0.01;% 0.01; % diffusion coefficient
        gamma = D ./ (1 -5*D); % anti-diffusion coefficient
        
        imFileName = 'circles.png'; %uigetfile(); % 'cameraman.tif';%'../figures/File0493.png';%
        % load image
        im = imread(imFileName);
        
        
        
    case 'longterm'
        
        invertColours = 0; % if true invert the colours of the image before starting
        nVisibleParticles = 0; % these are particles superposed to the image to
        % visualize the diffusion and anti-diffusion process: the image does not
        % interact with these particles, but they do interact with the image.
        
        onlyPlotAtTheEndOfCycle = 1;
        
        
        colourMapForParticles = jet(nVisibleParticles);
        imTitle = 'longterm_diff'; % title of saved figures
        imCounter = 1;
        figureFolder = 'figures';
        exportImageResolution = '-r300';
        isPlottingTrajectories = 1;
        clearTrajectoriesBeforeAntidiffusion = 0;
        
        iMax = 1; % the maximum value in the image at the beginning
        isBlurringWithLaplacian = 1; % using the Laplacian convolution is equivalent to taking a mean field approach
        isDeBlurringWithLaplacian = 1; % when deblurring, the laplacian convolution can produce images with negative values
        useMeanFieldForDeblurring = 1; % there isn't a useMeanFieldForBlurring because this would be equivalent to the Laplacian; in the case of deblurring
        % the further difference is whether the single particles that move
        % according to Weber's law do it probabilistically, or on average. If
        % isDeBlurringWithLaplacian, this useMeanField is ignored

%         reactionTime = 0; % if the reaction time is greater than zero, the
%         % de-blurring particles sense their environment at time t -
%         % reactionTime and move according to what they have sensed at time
%         % t
        avoidSaturationWhenDeblurring = 1; % if avoid saturation, then the particles won't move to a cell with more than iMax particles
        nBlurringSteps = 1; % 300
        nDeBlurringSteps = 1;
        nBlurringDeblurringCycles = 15000;
        
        D = 0.17;% 0.01; % diffusion coefficient
        gamma = 0.9;%0.9;% D ./ (1 -5*D); % anti-diffusion coefficient
        
        
        
        % % create an image for testing
        % im = ones(256)*50;
        if 0
            [X,Y] = meshgrid(-64:63);
            sigmaGaussianSeed = 10;
            baseline = 0
            im = (255 - baseline) * exp(- (X.^2 + Y.^2)/(sigmaGaussianSeed^2)) + baseline;
            % im(X.^2 + Y.^2 < 40^2) = 200;
        else
            im = rand(128) * 255;% rand(128)*128+127;
        end
        % im = 200*ones(size(X));
        % im(30,30) = 201;
        
end


if ~exist(figureFolder, 'dir') && isSavingFigures == 1
    mkdir(figureFolder);
end







allTrajectoriesX = [];
allTrajectoriesY = [];

if isPlottingTrajectories
    allTrajectoriesX = nan(nBlurringSteps + nDeBlurringSteps + 1, nVisibleParticles); % +1 is for the first image
    allTrajectoriesY = nan(nBlurringSteps + nDeBlurringSteps + 1, nVisibleParticles);
end

if iMax >= 10000 && useMeanFieldForDeblurring == 0 && isDeBlurringWithLaplacian == 0
    warndlg('This is going to be very slow because I use binornd for simulating particle movements. I suggest using a mean field approach');
end



selectedParticlePositions = nan(nVisibleParticles, 2);



if ndims(im) == 3
    imGray = (double(rgb2gray(im))*iMax/255); % TODO: I haven't really decided if I want to use double or int
else
    imGray = (double(im)*iMax/255);
end

% int seems appropriate for the binornd function below

if invertColours
    imGray = iMax - imGray;
end




% the size of the image on top of which the particles move is given by:
imWidth = size(imGray, 2);
imHeight = size(imGray,1);


% Get a certain number of particles for visualization
cumsumGrayLevels = reshape(cumsum(imGray(:)), imHeight, imWidth);

if isfinite(nVisibleParticles) && ceil(cumsumGrayLevels(end)) >= 1
    selectedParticles = randi(ceil(cumsumGrayLevels(end)), [nVisibleParticles, 1]);
    for jj = 1:nVisibleParticles
        [selectedParticlePositions(jj, 1),selectedParticlePositions(jj, 2)] = find(cumsumGrayLevels >= selectedParticles(jj), 1, 'first');
    end
end

% plot some particles on top of the image drawn from the distribution of
% pixel values

cla;
imagesc(imGray);% image(imGray*255/iMax)
colormap(gray)
caxis([0, iMax]);
hold on;
for pp = 1:nVisibleParticles
    plot(selectedParticlePositions(pp,2), selectedParticlePositions(pp,1), '.', 'MarkerSize', 22, 'Color', colourMapForParticles(pp,:));
end
axis off % not even necessary in this case
axis equal;
set(gca,'Units','normalized','Position',[0 0 1 1])
% print
% truesize(gcf,[imHeight imWidth]);


if isPlottingTrajectories && nVisibleParticles > 0
    allTrajectoriesX(imCounter,:) = selectedParticlePositions(:,2);
    allTrajectoriesY(imCounter,:) = selectedParticlePositions(:,1);
end

drawnow;



if isSavingFigures
    figureFileName = fullfile(figureFolder, sprintf('%s%05d.png',imTitle, imCounter));
    print(gcf, '-dpng', exportImageResolution, figureFileName);
    % set(gcf, 'InvertHardCopy', 'off');
    % print(gcf, '-depsc2', '-tiff', exportImageResolution, [figureFileName(1:end-3), 'eps']);
    % print(gcf, '-dtiff', exportImageResolution, [figureFileName(1:end-3), 'tiff']);
    
end



imCounter = imCounter + 1;


for ii = 1:nBlurringDeblurringCycles
    for jj = 1:nBlurringSteps
        imGray = perform_one_blur_step(imGray, D, isBlurringWithLaplacian);
        selectedParticlePositions = perform_diffusion_sel_particles(selectedParticlePositions, imWidth, imHeight, D);
        
        
%         sumImGrayBetweenDiffAndAntiDiff = (sum(imGray(:)))
        
        
        if onlyPlotAtTheEndOfCycle == 0
            
            
            cla;
            imagesc(imGray);% image(imGray*255/iMax)
            colormap(gray)
            caxis([0, iMax]);
            
            hold on;
            for pp = 1:nVisibleParticles
                plot(selectedParticlePositions(pp,2), selectedParticlePositions(pp,1), '.', 'MarkerSize', 22, 'Color', colourMapForParticles(pp,:));
            end
            axis off % not even necessary in this case
            axis equal;
            set(gca,'Units','normalized','Position',[0 0 1 1])
            % print
            % truesize(gcf,[imHeight imWidth]);
            
            if isPlottingTrajectories && nVisibleParticles > 0
                allTrajectoriesX(imCounter,:) = selectedParticlePositions(:,2);
                allTrajectoriesY(imCounter,:) = selectedParticlePositions(:,1);
                for pp = 1:nVisibleParticles
                    plot(allTrajectoriesX(:,pp), allTrajectoriesY(:,pp), '-', 'LineWidth', 2, 'Color', colourMapForParticles(pp,:));
                end
            end
            
            drawnow;
            
            if isSavingFigures
                figureFileName = fullfile(figureFolder, sprintf('%s%05d.png',imTitle, imCounter));
                print(gcf, '-dpng', exportImageResolution, figureFileName);
                % set(gcf, 'InvertHardCopy', 'off');
                % print(gcf, '-depsc2', '-tiff', exportImageResolution, [figureFileName(1:end-3), 'eps']);
                % print(gcf, '-dtiff', exportImageResolution, [figureFileName(1:end-3), 'tiff']);
                
            end
            
            
            
            imCounter = imCounter + 1;
        end
        
        

    end
    
    
    
    if isPlottingTrajectories && nVisibleParticles > 0
        if clearTrajectoriesBeforeAntidiffusion
            allTrajectoriesX = nan(nBlurringSteps + nDeBlurringSteps + 1, nVisibleParticles); % +1 is for the first image
            allTrajectoriesY = nan(nBlurringSteps + nDeBlurringSteps + 1, nVisibleParticles);
        end
    end
    
    
   
    
    if onlyPlotAtTheEndOfCycle == 0
        % Plot one more time, but without the trajectories
        cla;
        imagesc(imGray);% image(imGray*255/iMax)
        colormap(gray)
        caxis([0, iMax]);
        
        hold on;
        for pp = 1:nVisibleParticles
            plot(selectedParticlePositions(pp,2), selectedParticlePositions(pp,1), '.', 'MarkerSize', 22, 'Color', colourMapForParticles(pp,:));
        end
        axis off % not even necessary in this case
        axis equal;
        set(gca,'Units','normalized','Position',[0 0 1 1])
        % print
        % truesize(gcf,[imHeight imWidth]);
        
        if isPlottingTrajectories && nVisibleParticles > 0
            allTrajectoriesX(imCounter,:) = selectedParticlePositions(:,2);
            allTrajectoriesY(imCounter,:) = selectedParticlePositions(:,1);
            for pp = 1:nVisibleParticles
                plot(allTrajectoriesX(:,pp), allTrajectoriesY(:,pp), '-', 'LineWidth', 2, 'Color', colourMapForParticles(pp,:));
            end
        end
        
        drawnow;
        
        if isSavingFigures
            figureFileName = fullfile(figureFolder, sprintf('%s%05d.png',imTitle, imCounter));
            print(gcf, '-dpng', exportImageResolution, figureFileName);
            % set(gcf, 'InvertHardCopy', 'off');
            % print(gcf, '-depsc2', '-tiff', exportImageResolution, [figureFileName(1:end-3), 'eps']);
            % print(gcf, '-dtiff', exportImageResolution, [figureFileName(1:end-3), 'tiff']);
            
        end
        
        
        
        imCounter = imCounter + 1;
        
        
    end
    
    
    
    % pause(1);
    
    % apply Weber's law anti-diffusion
    for jj = 1:nDeBlurringSteps
        % imGray = perform_one_blur_step(imGray, D/2, isBlurringWithLaplacian);
        [imGray, pLeft, pRight, pUp, pDown] = perform_one_deblur_step(imGray, gamma, isDeBlurringWithLaplacian, useMeanFieldForDeblurring, avoidSaturationWhenDeblurring, iMax);
        selectedParticlePositions = perform_antidiffusion_sel_particles(selectedParticlePositions, imWidth, imHeight, pLeft, pRight, pUp, pDown);
        
        
        jj
        sumImGray = (sum(imGray(:)))
        minImGray = min(imGray(:))
        maxImGray = max(imGray(:))
        cla;
        imagesc(imGray);% image(imGray*255/iMax)
        colormap(gray)
        caxis([0, iMax]);
        
        hold on;
        for pp = 1:nVisibleParticles
            plot(selectedParticlePositions(pp,2), selectedParticlePositions(pp,1), '.', 'MarkerSize', 22, 'Color', colourMapForParticles(pp,:));
        end
        axis off % not even necessary in this case
        axis equal;
        set(gca,'Units','normalized','Position',[0 0 1 1])
        % print
        % truesize(gcf,[imHeight imWidth]);
        
        if isPlottingTrajectories && nVisibleParticles > 0
            allTrajectoriesX(imCounter,:) = selectedParticlePositions(:,2);
            allTrajectoriesY(imCounter,:) = selectedParticlePositions(:,1);
            for pp = 1:nVisibleParticles
                plot(allTrajectoriesX(:,pp), allTrajectoriesY(:,pp), '-', 'LineWidth', 2, 'Color', colourMapForParticles(pp,:));
            end
        end
        
        drawnow;
        
        if isSavingFigures
            % if mod(ii, 50) == 0
            figureFileName = fullfile(figureFolder, sprintf('%s%05d.png',imTitle, imCounter));
            print(gcf, '-dpng', exportImageResolution, figureFileName);
            % set(gcf, 'InvertHardCopy', 'off');
            % print(gcf, '-depsc2', '-tiff', exportImageResolution, [figureFileName(1:end-3), 'eps']);
            % print(gcf, '-dtiff', exportImageResolution, [figureFileName(1:end-3), 'tiff']);
            % end
        end
        
        
        
        imCounter = imCounter + 1;
    end
    
    % imGray = imGray .* (1+0.0001*randn(size(imGray))); % add some randomness
    
end

