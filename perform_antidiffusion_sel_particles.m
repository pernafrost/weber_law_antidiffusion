function selectedParticlePositions = perform_antidiffusion_sel_particles(selectedParticlePositions, imWidth, imHeight, pLeft, pRight, pUp, pDown)

% Notice that in this function I am not taking into account the
% avoidSaturationWhenDeblurring parameter!

nVisibleParticles = size(selectedParticlePositions,1);


% also move selected particles
% each visible particle climbs the gradient following the local weber
% contrast
indexSelectedParticlePositions = sub2ind(size(pLeft), selectedParticlePositions(:,1), selectedParticlePositions(:,2));
randomNumberSelectedParticles = rand(nVisibleParticles,1);

particleDirections = cumsum([pLeft(indexSelectedParticlePositions), pRight(indexSelectedParticlePositions), ...
    pUp(indexSelectedParticlePositions), pDown(indexSelectedParticlePositions)], 2);

% particles moving to the left are those that are moving and for which
% the random number, scaled to the total weberContrast, indicates left
selectedMovingLeft = find(randomNumberSelectedParticles <= particleDirections(:,1));
% particles moving to right are those that are moving and that are not
% moving to the left, and for which the Weber contrast indicates (left or right)
selectedMovingRight = setdiff(find(randomNumberSelectedParticles <= particleDirections(:,2)), selectedMovingLeft);
selectedMovingUp = setdiff(find(randomNumberSelectedParticles <= particleDirections(:,3)), union(selectedMovingLeft, selectedMovingRight));
selectedMovingDown = setdiff(find(randomNumberSelectedParticles <= particleDirections(:,4)), union(union(selectedMovingLeft, selectedMovingRight), selectedMovingUp));

if ~isempty(selectedMovingLeft)
    selectedParticlePositions(selectedMovingLeft,2) = mod(selectedParticlePositions(selectedMovingLeft,2) + imWidth - 2, imWidth) +1;
end
if ~isempty(selectedMovingRight)
    selectedParticlePositions(selectedMovingRight,2) = mod(selectedParticlePositions(selectedMovingRight,2), imWidth) +1;
end
if ~isempty(selectedMovingUp)
    selectedParticlePositions(selectedMovingUp,1) = mod(selectedParticlePositions(selectedMovingUp,1) + imHeight - 2, imHeight) +1;
end
if ~isempty(selectedMovingDown)
    selectedParticlePositions(selectedMovingDown,1) = mod(selectedParticlePositions(selectedMovingDown,1), imHeight) +1;
end