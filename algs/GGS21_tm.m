function [A_pr, ts, errCurve_all] = GGS21_tm(y, x, opts)
    % Gerchberg-Saxton 2-1 algorithm to retrieve complex TM. 
    % The forward model is y = abs(A*x), from y and x, it returns the estimate of A
	% Huang, G., et al. (2020). "Generalizing the Gerchberg–Saxton algorithm for retrieving complex optical transmission matrices." Photonics Research 9(1).
	
    saveErr = (nargout > 2);
    
    [m, p] = size(y);
    [n, ~] = size(x);
    
    %% Initialization
    y = y';
    iters = opts.iters;
    num_worker = opts.nWorker;
    useGPU = opts.useGPU;
    
    if ~isfield(opts, 'ratio'); opts.ratio = 0.89; end
    ratio = opts.ratio;

    if useGPU
        x = gpuArray(x');
    else
        x = x';
    end
    xinv = Tikinv(x); % xinv = x';
    
    Xfun  = @(t) x  * t;
    Xinvfun = @(t) xinv * t; 

    A_pr = (randn(n, m, 'single') + 1i*randn(n, m, 'single'))/sqrt(2);
    Y = y.*exp(1i*2*pi*rand(p, m, 'single')); % Phase initialization

    %% GGS2-1 algorithm 
    errCurve = zeros(m, iters, 'single');

    if useGPU
        nParallel = floor(m / num_worker); 
        ts = zeros(iters, ceil(m/num_worker));
    
        ii = 1;
        for i = 1:nParallel:m
            t0 = tic;
            iplus = min(i+nParallel-1, m);   
            Ysub = gpuArray(Y(:, i : iplus));
            ysub = gpuArray(y(:, i : iplus)); % target amplitude
            t=1;
            while t <= iters    
                curr_Asub = Xinvfun(Ysub);
                Ysub_iter = Xfun(curr_Asub);
                if t <= round(ratio * iters) % GS-2
                    Ysub = ysub.^2.*exp(1i*angle(Ysub_iter)); 
                else % GS-1         
                    Ysub = ysub.*exp(1i*angle(Ysub_iter));
                end
                if saveErr; errCurve(i : iplus, t) = gather(vecnorm(abs(Ysub_iter) - ysub, 2, 1)); end
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
        while t <= iters    
            A_pr = Xinvfun(Y);
            Y_iter = Xfun(A_pr);
            if t <= round(ratio * iters) % GS-2
                Y = y.^2.*exp(1i*angle(Y_iter)); 
            else % GS-1         
                Y = y.*exp(1i*angle(Y_iter));
            end
            if saveErr; errCurve(:, t) = vecnorm(y-abs(Y_iter), 2, 1); end
            ts(t) = toc(t0);
            t = t+1;
        end
    end

    if saveErr; errCurve_all = mean(errCurve, 1)'; end
    A_pr = A_pr';
    
    % normalization to std=1
    infInd = find(abs(A_pr)>1e10);
    A_pr(infInd) = exp(1i*angle(A_pr(infInd)));
    A_pr = 1/std(A_pr, 1, 'all') * A_pr;
    
end