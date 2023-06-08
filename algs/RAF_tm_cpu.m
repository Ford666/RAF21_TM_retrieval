function [A_pr, ts, errCurve, corrCurve] = RAF_tm_cpu(y, x, opts)
    % Reweighted amplitude flow model to retrieve TM. 
    % The forward model is Y = abs(A * X).^2, from Y and X, it returns the estimate of A
    % CPU version, with either iteration number or running time as the stop criterion
    
    saveErr = (nargout > 2);
    saveCorr = (nargout > 3);
    
    [m, p] = size(y);
    [n, ~] = size(x);
    
    %% Initialization
    npower_iter = 30;                                                                               % Number of power iterations 
    beta = 5;
    gamma = 0.5; 
    mu = opts.mu;

    y = y'; ymag = sqrt(y);

    x = x';
    xH = x';
    
    Xfun  = @(I) x  * I;
    Xinvfun = @(Y) xH * Y;
   
    
    A_pr = (randn(n, m, 'single') + 1i*randn(n, m, 'single'))/sqrt(2); 
    A_pr = A_pr./vecnorm(A_pr, 2, 1);                                                               % initial guess
    
    normest = sqrt(mean(y, 1));                                                                     % norm to scale eigenvector
    ysort = sort(y, 1, 'ascend');                                                                   % sorting at ascending order
    ythresh = ysort(round(p / 1.3), :);
    ind = y >= ythresh;
    
    %% Weighted max. correlation initialization
    weights = (ymag .* ind).^gamma;                                                                 % remove and reweight the components
    for tt = 1:npower_iter                                                                          % power iteration     
        A_pr = Xinvfun( weights.* (Xfun(A_pr)) ); 
        A_pr = A_pr./vecnorm(A_pr, 2, 1);
    end
    
    A_pr = normest .* A_pr;                                                                         % Apply scaling
    
    
    %% Gradient-descent loop 
    if isfield(opts, 'iters') 
        iters = opts.iters;
        errCurve = zeros(iters, 1);
        corrCurve = zeros(iters, 1);
        ts = zeros(iters, 1);

        t0 = tic;
        compt = 1;
        
        % RAF
        while compt <= iters   
            curr_y = Xfun(A_pr);
            weight = abs(curr_y) ./ (abs(curr_y)+beta*ymag); % gradient weight
            grad = Xinvfun( (curr_y - ymag .* curr_y./abs(curr_y)) .* weight) / p;
            A_pr = A_pr - mu * grad;
            if saveErr; errCurve(compt) = sum(vecnorm(ymag - abs(curr_y), 2, 1)); end  % least-square measurement error
            if saveCorr; corrCurve(compt) = sum(vecCorr(ymag, abs(curr_y), 1)); end    % correlation coefficient to the true measurement
            ts(compt) = toc(t0);
            compt = compt+1;
        end

    elseif isfield(opts, 'times') 
        errCurve = [];
        corrCurve = [];
        ts = [0];

        t0 = tic;

        % RAF
        while ts(end) <= opts.times    
            curr_y = Xfun(A_pr);
            weight = abs(curr_y) ./ (abs(curr_y)+beta*ymag); % gradient weight
            grad = Xinvfun( (curr_y - ymag .* curr_y./abs(curr_y)) .* weight) / p;
            A_pr = A_pr - mu * grad;
            if saveErr; errCurve = [errCurve, sum(vecnorm(ymag - abs(curr_y), 2, 1))]; end  % least-square measurement error
            if saveCorr; corrCurve = [corrCurve, sum(vecCorr(ymag, abs(curr_y), 1))]; end    % correlation coefficient to the true measurement
            ts = [ts, toc(t0)];

        end

    end

    if saveErr; errCurve = errCurve'./m; end
    if saveCorr; corrCurve = corrCurve'./m; end
    A_pr = A_pr';

    % normalization to std=1
    infInd = find(abs(A_pr)>1e10);
    A_pr(infInd) = exp(1i*angle(A_pr(infInd)));
    A_pr = 1/std(A_pr, 1, 'all') * A_pr;
end