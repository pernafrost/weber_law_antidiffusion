function selectedParticlePositions = perform_diffusion_sel_particles(selectedParticlePositions, imWidth, imHeight, D)


nVisibleParticles = size(selectedParticlePositions,1);

% simulate the movement of selected particles
% each visible particle moves randomly with probability D*4, and if it
% moves it picks randomly one direction out of four possible directions
movingDirectionSelectedParticles = (rand(nVisibleParticles, 1) <= D*4) .* randi(4, [nVisibleParticles,1]);
selectedMovingLeft = find(movingDirectionSelectedParticles == 1);
selectedMovingRight = find(movingDirectionSelectedParticles == 2);
selectedMovingUp = find(movingDirectionSelectedParticles == 3);
selectedMovingDown = find(movingDirectionSelectedParticles == 4);

selectedParticlePositions(selectedMovingLeft,2) = mod(selectedParticlePositions(selectedMovingLeft,2) + imWidth - 2, imWidth) +1;
selectedParticlePositions(selectedMovingRight,2) = mod(selectedParticlePositions(selectedMovingRight,2), imWidth) +1;
selectedParticlePositions(selectedMovingUp,1) = mod(selectedParticlePositions(selectedMovingUp,1) + imHeight - 2, imHeight) +1;
selectedParticlePositions(selectedMovingDown,1) = mod(selectedParticlePositions(selectedMovingDown,1), imHeight) +1;

