function [z, ts, errCurve, corrCurve] = RAF21_objRe(y, A, opts)
    % RAF21 model to reconstruct object from speckle pattern measurements. 
    % The forward model is y = abs(A * x), from y and A, it returns the estimate of x
    % RAF2-1, adopted from RAF
    
    saveErr = (nargout > 2) && (isfield(opts, 'obj'));
    saveCorr = (nargout > 3) && (isfield(opts, 'obj'));

    t0 = tic;
    [m, n] = size(A);
    
    Ainv = A';

    ymag = sqrt(y); 
       
    Afun  = @(t) A * t;
    Ainvfun = @(t) Ainv * t;  
    
    %% Initialization
    npower_iter = 100; % Number of power iterations 
    beta = 10; %real: 10 / complex: 5
    gamma = 0.5;
    mu = 2.2;
    iters = opts.iters;
    RAF21_flag = 0;
    if isfield(opts, 'ratio'); ratio = opts.ratio; RAF21_flag = 1; else; ratio=2/3; end

    if strcmp(opts.signalType, 'complex') 
        z0 = (randn(n, 1, 'single') + 1i*randn(n, 1, 'single'))/sqrt(2);
    else
        z0 = randn(n, 1, 'single');
    end
    z0 = z0./norm(z0); %initial guess
    
    normest = sqrt(mean(y, 1)); % norm to scale eigenvector
%     normest = (abs(Afun(z0))'*ymag) / (abs(Afun(z0))'*abs(Afun(z0)));  % for non-Gaussian signal

    ysort = sort(y, 1, 'ascend');  % sorting at ascending order
    ythresh = ysort(round(m / 1.3), :);
    ind = y >= ythresh;
    
    %% Weighted max. correlation initialization
    weight = (ymag .* ind).^gamma;

    for tt = 1:npower_iter                % power iteration     
        z0 = Ainvfun( weight.* Afun(z0) ); 
        z0 = z0./norm(z0);
    end
    
    z0 = normest .* z0;  % Apply scaling
    
    
    %% Gradient-descent loop 
    if saveErr; errCurve = zeros(iters, 1, 'single'); end
    if saveCorr; corrCurve = zeros(iters, 1, 'single'); end
    if isfield(opts, 'obj'); x = reshape(opts.obj, [], 1); end

    z = z0;
    t = 1;
    while t <= iters    
        yz = Afun(z);
        weight = 1 ./ (1+beta./(abs(yz) ./ ymag)); % gradient weight
        if t <= ceil(ratio*iters) && RAF21_flag
            grad = 1/m * Ainvfun( (yz - ymag.^2 .* yz./abs(yz)) .* weight);  
        else
            grad = 1/m * Ainvfun( (yz - ymag .* yz./abs(yz)) .* weight);  
        end
        z = z - mu * grad;
        if saveErr; errCurve(t, :) = norm(ymag - abs(yz), 'fro') / norm(ymag, 'fro'); end
        if saveCorr; corrCurve(t, :) = corr(real(x), real(exp(-1i * angle(trace(x' * z))) * z)); end
        t=t+1;
    end
    if strcmp(opts.signalType, 'complex') 
        z = angle(z);
    end
    
    ts = toc(t0);
end