%
% %% here is a test:
% %singleOne = zeros(3); singleOne(2,2) = 1;
% singleOne = zeros(5); singleOne(3,3) = 1;
% D = 1/6;
% gamma = 0.3
% lDiff = [0, 1, 0; 1, -4, 1; 0, 1, 0]; % this is one version of discrete laplace operator
% singleOneDiff = singleOne + convn(singleOne, lDiff, 'same')*D % diffusion
% lF = convn(singleOneDiff, - lDiff, 'same')*(1/D) % gradient climbing; this is the filter involving laplacian and bilaplacian
% b = singleOneDiff - gamma * D * lF % I don't remember why I was testing this. Probably it is to see if it produces branching starting from an isotropic distribution?
% b = singleOne + convn(singleOne, lDiff, 'same')*D - gamma * D * lF % this is the same as above, but involves the simultaneous diffusion and contrast climbing
%
%
%
% % sum of b along the vertical and horizontal direction:
% sum(b((end +1)/2,:))
% sum(sum(b.*eye(size(b))))
%


%% Here is the actual code

clear all;
close all;

initialCondition = 'rand'; %'mixed'; % 'gauss'; %


rng(10000); % With this line I use the same generator


isSavingFigures = 1;
simulatePattern = 1;

figureFolder = 'figures_long_term';

if ~exist(figureFolder, 'dir') && isSavingFigures == 1
    mkdir(figureFolder);
end

w = 128; % side of the lattice

% boundaryConditions = 'periodic';
b = 1; % size of periodic boundary
wp = w + 2*b; % add periodic boundary at the beginning and at the end of each lattice dimension
wbeg = b+1; % beginning of actual density landscape (not considering periodic boundaries)
wend = w+b; % end of actual density landscape (not considering periodic boundaries)

% apply gradient response
s = 0.000; % this controls an additional increment of density
singleOne = zeros(5); singleOne(3,3) = 1;
lDiff = [0, 1, 0; 1, -4, 1; 0, 1, 0]; % this is one version of discrete laplace operator



% lF = convn(lDiff, lDiff) + 6 * convn(lDiff, singleOne)

% Up to gamma =1.5 D =0.1
% from gamma 0.35 D 0.05
% best gamma 0.9 D 0.1666
% gamma = 0.5 D = 0.1

for gamma = 0.9; %0.1:0.1:1; % antidiffusion coefficient
    for D = 0.13:0.04:0.17;% 0.05:0.05:0.4; % diffusion coefficient
        
        % this tests a single large vector
        singleOneDiff = singleOne + convn(singleOne, lDiff, 'same')*D; % diffusion
        % lF = convn(singleOneDiff, - lDiff, 'same')*(1/D); % gradient climbing; this is the filter involving laplacian and bilaplacian
        % f = convn(singleOne, lDiff, 'same')*D - gamma * D * lF;
        
        % sinkCoordinates = [wend, wend]; % these are coordinates of a sink that always remains empty
        
        
        % initialize lattice
        x=zeros([wp, wp]); % the size of the lattice is its real size plus the periodic boundaries
        
        switch initialCondition
            case 'rand'
                % either put random noise everywhere
                x(wbeg:wend, wbeg:wend) = rand([w, w]);
                
            case 'gauss'
                % or initiate with a gaussian somewhere
                x(wbeg+w/2, wbeg+w/2) = 1; % add a single germ somewhere
                for jj = 1:200
                    x = x + D*convn(x,lDiff,'same');
                    x = wextend('2D','ppd',x(wbeg:wend,wbeg:wend),[b, b]);
                end
                xmax = max(x(:));
                x = x/xmax;
            case 'mixed'
                x(wbeg+w/2, wbeg+w/2) = 1; % add a single germ somewhere
                for jj = 1:200
                    x = x + D*convn(x,lDiff,'same');
                    x = wextend('2D','ppd',x(wbeg:wend,wbeg:wend),[b, b]);
                end
                xmax = max(x(:));
                x = x/xmax*0.99;
                x(wbeg:wend, wbeg:wend) = x(wbeg:wend, wbeg:wend) + rand([w, w])*0.01;

        end
        
        % switch boundaryConditions
        %     case 'periodic'
        % x(wend+1:end,:) = x(wbeg:wbeg-1+b,:); % the actual vector is 2:100 and I apply "periodic boundary conditions
        % x(1:b) = x(wend -b+1:wend);
        x = wextend('2D','ppd',x(wbeg:wend,wbeg:wend),[b, b]);
        %     case 'dirichlet'
        %         x(wend+1:end) = T;% x(wend); % the actual vector is 2:100 and I apply "periodic boundary conditions
        %         x(1:b) = T;% x(wbeg);
        %     otherwise
        %         errordlg('Boundary conditions can be ''periodic'' or ''dirichlet''', 'Error in the code');
        %         return;
        % end
        
        figure, imagesc(x(wbeg:wend, wbeg:wend)); colormap(gray);
        
        lF = convn(singleOneDiff, - lDiff, 'same')*(1/D) % gradient climbing; this is the filter involving laplacian and bilaplacian
        bMatrix = singleOne + convn(singleOne, lDiff, 'same')*D - gamma * D * lF % this is the same as above, but involves the simultaneous diffusion and contrast climbing
        
        if simulatePattern
            for jj = 1:15000
                
                % x(sinkCoordinates) = 0; % the sink is always at zero
                % compute the skeleton of empty space
                
                % skeleton = bwmorph(x<0.5, 'skel', 'Inf');
                % skeletonDiffused = 1 - (skeleton + D*convn(skeleton,lDiff,'same'));
                
                % x = x - gamma * skeletonDiffused .* (x.*(1-x)) .* convn(x,lF, 'same') + D*convn(x,lDiff,'same') + x.*(1-x).*rand(size(x))*s .* skeletonDiffused;
                x1 = x + D*convn(x,lDiff,'same');
                x1 = wextend('2D','ppd',x1(wbeg:wend,wbeg:wend),[b, b]);
                x = x1 - gamma * (x1.*(1-x1)) .* convn(x1, lDiff, 'same') + x1.*(1-x1).*rand(size(x1))*s;
                
                
                sum(x(:)) % sum should remain constant if s = 0
                
                % x = x - (x.*(1-x)) .* convn(x, f, 'same');
                
                
                % There are some negative numbers appearing. For the moment I set them to
                % zero
                % x = max(x, 0);
                
                
                % convn(x,lF, 'same')
                % x(wbeg:wend,wbeg:wend).*(1-x(wbeg:wend,wbeg:wend))
                
                x = wextend('2D','ppd',x(wbeg:wend,wbeg:wend),[b, b]);
                
                if mod(jj, 100) == 0
                    
%                     h=pcolor(x(wbeg:wend, wbeg:wend));
%                     set(h,'edgecolor','none');
%                     shading flat;
                    imagesc(x(wbeg:wend, wbeg:wend));
                    axis square; axis tight;
                    colormap(gray);
                    caxis([0, 1]);
                    axis off;
                    
                    
                    
                    % fprintf([repmat('%12.8f',1,size(x(wbeg:wend,wbeg:wend),2)) '\n'],x(wbeg:wend, wbeg:wend))
                    drawnow; pause(0.1);
                    jj
                end
                
                
            end
            
            
            maxX = max(x(:))
            minX = min(x(:))
            
            % title(sprintf('gamma: %0.4f; D: %0.4f', gamma, D));
            
            if isSavingFigures
                figureFileName = fullfile(figureFolder, sprintf('sim2D_time_sep_gamma%dover1000_D%dover1000.png',round(gamma*1000), round(D*1000)));
                print(gcf, '-dpng', '-r300', figureFileName);
                print(gcf, '-depsc2', '-tiff', '-r300', [figureFileName(1:end-3), 'eps']);
                imwrite(x(wbeg:wend, wbeg:wend), [figureFileName(1:end-3), 'tif']);
                
                % also visualise and save Fourier amplitude
                x1 = x(wbeg:wend, wbeg:wend);
                x1F = (fft2(x1));
                x1F(1,1)=0; % remove DC
                x1F = fftshift(x1F);
                figure, h = pcolor(-64.5:62.5, -64.5:62.5, abs(x1F)); 
                set(gca,'FontSize', 34, 'TickDir', 'out'); box off;
                colormap(gray); shading flat; set(h,'edgecolor','none');
                axis equal; axis tight;
                xlabel('Freq. (cycles / image)');
                ylabel('Freq. (cycles / image)');
                
                figureFileName2 = [figureFileName(1:end-4), '_amp_spectrum.png']
                print(gcf, '-dpng', '-r300', figureFileName2);
                print(gcf, '-depsc2', '-tiff', '-r300', [figureFileName2(1:end-3), 'eps']);
                imwrite(x1F, [figureFileName2(1:end-3), 'tif']);

            end
        end
        
        close all;
        % pause(60);
    end
end



