function [A_pr, ts, errCurve, corrCurve] = RAF_tm(y, x, opts)
    % Reweighted amplitude flow model to retrieve TM. 
    % The forward model is Y = abs(A * X).^2, from Y and X, it returns the estimate of A
    
    saveErr = (nargout > 2);
    saveCorr = (nargout > 3);
    
    [m, p] = size(y);
    [n, ~] = size(x);
    
    %% Initialization
    npower_iter = 30;                                                                               % Number of power iterations 
    beta = 5;
    gamma = 0.5; 
    mu = opts.mu;
    iters = opts.iters;
    num_worker = opts.nWorker;

    if ~isfield(opts, 'useGPU'); opts.useGPU = 0; end
    if ~isfield(opts, 'ratio'); opts.ratio = 2/3; end
    useGPU = opts.useGPU;
    ratio = opts.ratio;

    y = y'; ymag = sqrt(y);
    if useGPU
        x = gpuArray(x');
    else
        x = x';
    end
    xH = x';
    
    Xfun  = @(I) x  * I;
    Xinvfun = @(Y) xH * Y;
   
    
    if useGPU
        A_pr = (randn(n, m, 'single', 'gpuArray') + 1i*randn(n, m, 'single', 'gpuArray'))/sqrt(2); 
    else
        A_pr = (randn(n, m, 'single') + 1i*randn(n, m, 'single'))/sqrt(2); 
    end
    A_pr = A_pr./vecnorm(A_pr, 2, 1);                                                               % initial guess
    
    normest = sqrt(mean(y, 1));                                                                     % norm to scale eigenvector
    if useGPU; normest = gpuArray(normest); end
    ysort = sort(y, 1, 'ascend');                                                                   % sorting at ascending order
    ythresh = ysort(round(p / 1.3), :);
    ind = y >= ythresh;
    
    %% Weighted max. correlation initialization
    weights = (ymag .* ind).^gamma;                                                                 % remove and reweight the components
    if useGPU; weights = gpuArray(weights); end
    for tt = 1:npower_iter                                                                          % power iteration     
        A_pr = Xinvfun( weights.* (Xfun(A_pr)) ); 
        A_pr = A_pr./vecnorm(A_pr, 2, 1);
    end
    
    A_pr = gather(normest .* A_pr);                                                                         % Apply scaling
    
    
    %% Gradient-descent loop 
    errCurve = zeros(m, iters, 'single');
    corrCurve = zeros(m, iters, 'single');
                                                                                                   % parameter for RAF, real: 10 / complex: 5
    if useGPU
        nParallel = floor(m / num_worker);
        ts = zeros(iters, ceil(m/num_worker));

        ii = 1;
        for i = 1:nParallel:m
            t0 = tic;
            iplus = min(i+nParallel-1, m);   
            ysub = gpuArray(ymag(:, i : iplus)+1.0e-30);
            curr_Asub = gpuArray(A_pr(:, i : iplus));
            t = 1;

            % RAF2-1
            while t <= iters   
                curr_y = Xfun(curr_Asub);
                weight = abs(curr_y) ./ (abs(curr_y)+beta*ysub); % gradient weight
                grad = Xinvfun( (curr_y - ysub .* curr_y./abs(curr_y)) .* weight) / p;
                curr_Asub = curr_Asub - mu * grad;
                if saveErr; errCurve(i : iplus, t) = gather(vecnorm(ysub - abs(curr_y), 2, 1)).'; end  % least-square measurement error
                if saveCorr; corrCurve(i : iplus, t) = gather(vecCorr(ysub, abs(curr_y), 1)).'; end    % correlation coefficient to the true measurement
                ts(t, ii) = toc(t0);
                t = t+1;
            end  
            A_pr(:, i : iplus) = gather(curr_Asub);  
            ii = ii+1;
        end
        ts = sum(ts, 2);
    else
        ts = zeros(iters, 1);
        t0 = tic;
        t = 1;

        % RAF2-1
        while t <= iters   
            curr_y = Xfun(A_pr);
            weight = abs(curr_y) ./ (abs(curr_y)+beta*ymag); % gradient weight
            grad = Xinvfun( (curr_y - ymag .* curr_y./abs(curr_y)) .* weight) / p;
            A_pr = A_pr - mu * grad;
            if saveErr; errCurve(:, t) = vecnorm(ymag - abs(curr_y), 2, 1).'; end  % least-square measurement error
            if saveCorr; corrCurve(:, t) = vecCorr(ymag, abs(curr_y), 1).'; end    % correlation coefficient to the true measurement
            ts(t) = toc(t0);
            t = t+1;
        end     
    end

    if saveErr; errCurve = mean(errCurve, 1)'; end
    if saveCorr; corrCurve = mean(corrCurve, 1)'; end
    A_pr = A_pr';

    % normalization to std=1
    infInd = find(abs(A_pr)>1e10);
    A_pr(infInd) = exp(1i*angle(A_pr(infInd)));
    A_pr = 1/std(A_pr, 1, 'all') * A_pr;

end