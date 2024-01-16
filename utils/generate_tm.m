function TM = generate_tm(M, N)
    %generate_tm generates a Transmission Matrix (TM) that relates one 
    %unit input pixel to its resultant speckle field
    % Parameters:
    %               P: height of Transmission Matric
    %               N: weight of Transmission Matric
    % Returns:
    %               Transmission Matric

    % Adapted from https://github.com/comediaLKB/NMF_PR/blob/master/generate_tm.m
    
    % Pupil definition    
    sz_grains = 4;   %size of a single speckle grain
    if mod(round(sqrt(M)/sz_grains), 2)==0; n_grains  = round(sqrt(M)/sz_grains); else; n_grains  = round(sqrt(M)/sz_grains)+1; end%number of speckle grains
    
    [Y, X] = meshgrid(1:n_grains, 1:n_grains);
    pupil = (X-n_grains/2).^2 + (Y-n_grains/2).^2 < (n_grains/2)^2;
    pupil = single(pupil);
    
    % A speckle is modeled by a random phase in the Fourier space
    bruit = exp(2*1i*pi * rand(n_grains, n_grains, N, 'single'));  %spectral speckle 
    bruit = bruit .* pupil;              %imaging system CTF/pupil func
    
    % Fourier transform to go in the object space with zero padding for smoothing
    TM = zeros(M, N, 'single');

    for j = 1:N
        temp = fft2(fftshift(padarray(bruit(:, :, j), ...
            [sqrt(M)/2 - n_grains/2, sqrt(M)/2 - n_grains/2])));
        TM(:, j) = temp(:);
    end
    
end


