function [A_pr, ts, errCurve] = GGS21_tm_cpu(y, x, opts)
    % Gerchberg-Saxton 2-1 algorithm to retrieve complex TM. 
    % The forward model is y = abs(A*x), from y and x, it returns the estimate of A
	% Huang, G., et al. (2020). "Generalizing the Gerchberg–Saxton algorithm for retrieving complex optical transmission matrices." Photonics Research 9(1).
	% CPU version, with either iteration number or running time as the stop criterion

    saveErr = (nargout > 2);
    
    [m, p] = size(y);
    [n, ~] = size(x);
    
    %% Initialization
    y = y';

    if ~isfield(opts, 'ratio'); opts.ratio = 0.89; end
    ratio = opts.ratio;


    x = x';

    xinv = Tikinv(x); % xinv = x';
    
    Xfun  = @(t) x  * t;
    Xinvfun = @(t) xinv * t; 

    A_pr = (randn(n, m, 'single') + 1i*randn(n, m, 'single'))/sqrt(2);
    Y = y.*exp(1i*2*pi*rand(p, m, 'single')); % Phase initialization

    %% GGS2-1 algorithm 

    if isfield(opts, 'iters') 
        iters = opts.iters;
        errCurve = zeros(iters, 1);
        ts = zeros(iters, 1);
  
        t0 = tic;
        compt = 1;
        while compt <= iters    
            A_pr = Xinvfun(Y);
            Y_iter = Xfun(A_pr);
            if compt <= round(ratio * iters) % GS-2
                Y = y.^2.*exp(1i*angle(Y_iter)); 
            else % GS-1         
                Y = y.*exp(1i*angle(Y_iter));
            end
            if saveErr; errCurve(compt) = sum(vecnorm(y-abs(Y_iter), 2, 1)); end
            ts(compt) = toc(t0);
            compt = compt+1;
        end

    elseif isfield(opts, 'times') 
        errCurve = [];
        ts = [0];

        t0 = tic;

        while ts(end) <= opts.times    
            A_pr = Xinvfun(Y);
            Y_iter = Xfun(A_pr);
            if ts(end) <= round(ratio * opts.times) % GS-2
                Y = y.^2.*exp(1i*angle(Y_iter)); 
            else % GS-1         
                Y = y.*exp(1i*angle(Y_iter));
            end
            if saveErr; errCurve = [errCurve, sum(vecnorm(y-abs(Y_iter), 2, 1))]; end
            
            ts = [ts, toc(t0)];
        end

    end

    if saveErr; errCurve = errCurve'./m; end
    A_pr = A_pr';
    
    % normalization to std=1
    infInd = find(abs(A_pr)>1e10);
    A_pr(infInd) = exp(1i*angle(A_pr(infInd)));
    A_pr = 1/std(A_pr, 1, 'all') * A_pr;
    
end