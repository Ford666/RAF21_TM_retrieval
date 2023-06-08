function [A_pr, ts, errCurve, corrCurve] = RAF_tm_cpu_gradAccel(y, x, opts)
    % Reweighted amplitude flow model to retrieve TM. 
    % The forward model is Y = abs(A * X).^2, from Y and X, it returns the estimate of A
    % Wang, G., et al. (2018). "Phase retrieval via reweighted amplitude flow." IEEE Transactions on Signal Processing 66(11): 2818-2833.
    
    % Date: April 2023
    % Author: Shengfu Cheng
    % This version: running on CPU with time counting each iteration

    % 1) Adaptive stepsize via line search (Barzilai-Borwein method)
    % 2) Advanced gradient descent method (conjugate gradient, L-BFGS etc.)


    saveErr = (nargout > 2);
    saveCorr = (nargout > 3);
    
    muRAF = opts.mu;
    searchMethod = opts.searchMethod;
    if strcmpi(opts.searchMethod, 'NCG')
        if ~isfield(opts, 'beta')
            opts.beta = 'FR';
        end
        CGresetPeriod = 2;
    elseif strcmpi(opts.searchMethod, 'LBFGS')
        if ~isfield(opts, 'stores')
            opts.stores = 5;
        end
        BFGSresetPeriod = 20;
    end
    if ~isfield(opts, 'updateMu')
        opts.updateMu = 'False';
    end
    
    [m, p] = size(y);
    [n, ~] = size(x);
    
    y = y'; ymag = sqrt(y); 
    x = x';
    xH = x';
    
    Xfun  = @(I) x  * I;
    Xinvfun = @(Y) xH * Y;
    
    %% Initialization
    npower_iter = 30; % Number of power iterations 
    alp = 5;   %real: 10 / complex: 5
    ga = 0.5; 
    
    A_pr = (randn(n, m, 'single') + 1i*randn(n, m, 'single'))/sqrt(2); 
    A_pr = A_pr./vecnorm(A_pr, 2, 1); %initial guess
    
    normest = sqrt(mean(y, 1)); % norm to scale eigenvector
    ysort = sort(y, 1, 'ascend');  % sorting at ascending order
    ythresh = ysort(round(p / 1.3), :);
    ind = y >= ythresh;
    
    %% Weighted max. correlation initialization
    weights = (ymag .* ind).^ga;  %remove and reweight the components
    for tt = 1:npower_iter                     %power iteration     
        A_pr = Xinvfun( weights.* (Xfun(A_pr)) ); 
        A_pr = A_pr./vecnorm(A_pr, 2, 1);
    end
 
    A_pr = normest .* A_pr;  % Apply scaling
    
    %% Gradient-descent loop 
    if isfield(opts, 'iters') 
        iter_RAF = opts.iters;
        errCurve = zeros(iter_RAF, 1);
        corrCurve = zeros(iter_RAF, 1');
        
        ts = zeros(iter_RAF, 1);

        t0 = tic;
        compt = 1; compt0 = 1;
        
        while compt <= iter_RAF      
            curr_y = Xfun(A_pr);
            weight = abs(curr_y) ./ (abs(curr_y)+alp*ymag); % gradient weight
            grad = Xinvfun((curr_y - ymag.*curr_y./abs(curr_y)).*weight)/p;
            
            % update via search direction & stepsize
            switch searchMethod
                case 'SD'
                    searchDir = -grad;
                case 'NCG'
                    if (compt==1 || compt-compt0 == CGresetPeriod) 
                        searchDir0 = zeros(size(grad), 'single');
                        beta = zeros(1, m, 'single');
                        compt0 = compt;
                    else       
                        if strcmpi(opts.beta, 'HS')
                            Dg = grad - grad0;
                            beta = -real(sum(conj(grad).*Dg, 1)) ./ real(sum(conj(searchDir0).*Dg, 1));
                        elseif strcmpi(opts.beta, 'FR')
                            beta = sum(grad.^2, 1) ./ sum(grad0.^2, 1);
                        end
                    end
                    searchDir = -grad + beta .*searchDir0;

                case 'LBFGS'
                    searchDir = -grad;
                    iterSub = min(compt-compt0, opts.stores);

                    if (compt==1 || compt-compt0 == BFGSresetPeriod) 
                        % Perform LBFGS initialization
                        yVals = zeros(n, m, opts.stores, 'single');
                        sVals = zeros(n, m, opts.stores, 'single');
                        rhoVals = zeros(1, m, opts.stores, 'single');
                        compt0 = compt;
                    else
                        % Update LBFGS stored vectors
                        Dg = grad - grad0; 
                        sVals = cat(3, Dx, sVals(:, :, 1:opts.stores-1));
                        yVals = cat(3, Dg, yVals(:, :, 1:opts.stores-1));
                        rhoVals = cat(3, 1./real(sum(conj(Dg).*Dx, 1)), rhoVals(:, :, 1:opts.stores-1));
                    end  
                    
                    if iterSub > 0  % two loops to determine the search direction
                        alphas = zeros(1, m, iterSub, 'single');
                
                        % First loop
                        for j = 1 : iterSub
                            alphas(:,:,j) = rhoVals(:,:,j) .* real(sum(conj(sVals(:,:,j)) .* searchDir, 1));
                            searchDir = searchDir - alphas(:,:,j) .* yVals(:,:,j);
                        end
                        
                        % Scaling of search direction
                        gamma = real(sum(conj(Dg) .* Dx, 1)) ./ sum(conj(Dg) .* Dg, 1);
                        searchDir = gamma .* searchDir;
                        
                        % Second loop
                        for j = iterSub : -1 : 1
                            beta = rhoVals(:,:,j) .* real(sum(conj(yVals(:,:,j)) .* searchDir, 1));
                            searchDir = searchDir + (alphas(:,:,j) - beta) .* sVals(:,:,j);
                        end
                        
                        searchDir = 1./gamma .* searchDir;
                        searchDir = vecnorm(grad, 2, 1) ./ vecnorm(searchDir, 2, 1) .* searchDir;
                    end
            end
            if strcmp(opts.updateMu, 'True') || strcmp(searchMethod, 'NCG') || strcmp(searchMethod, 'LBFGS')
                grad0 = grad;
            end
            if strcmp(opts.updateMu, 'True') || strcmp(searchMethod, 'NCG')
                searchDir0 = searchDir;
            end
 
            % update the stepsize with Barzilai-Borwein adaptive method
            if strcmp(opts.updateMu, 'True')
                if mod(compt-1, 3)==0                           % reinitialization of stepsize
                    mu1 = muRAF;
                else                                           % update the stepsize
                    Ds = searchDir - searchDir0;
                    dotprod = real(sum(conj(Dx).*Ds, 1));
                    muS = vecnorm(Dx,2,1).^2 ./ dotprod;       % First BB stepsize rule
                    muM = dotprod ./ vecnorm(Ds,2,1).^2;       % Alternate BB stepsize rule
                    muM = max(muM, 0);
                    ind = find(2*muM <= muS);
                    mu1 = muM; mu1(ind) = muS(ind)-muM(ind)./2; 
                    ind = union(find(mu1<=0), find(isinf(mu1))); 
                    ind = union(ind, find(isnan(mu1)));
                    if numel(mu0) > 1
                        mu1(ind) = 1.5*mu0(ind);              % grows first, then shinks with backtracking
                    else
                        mu1(ind) = 1.5*mu0;
                    end
                end
                f0 = 0.5/p * vecnorm(abs(curr_y)-ymag, 2, 1).^2;
            else
                mu1 = muRAF;
            end
            

            % update the TM
            A_pr = A_pr + mu1 .* searchDir;
            curr_y = Xfun(A_pr);
            
            if strcmp(opts.updateMu, 'True')
                mu0 = mu1; 
            end
 
            % Backtracking to find the stepsize with Armijo-Goldstein condition
            if strcmp(opts.updateMu, 'True') && mod(compt-1, 3)~=0
                ind = 1:m;
                f1 = 0.5/p * vecnorm(abs(curr_y(:, ind))-ymag(:, ind), 2, 1).^2;
                backtrackCount = 0;
                while backtrackCount <= 20 
                    tmp = f0(:, ind) + 0.1 * mu0(:, ind) .* real(sum(conj(searchDir0(:, ind)).*grad0(:, ind), 1));  
                    ind = ind(f1>tmp);
                    if isempty(ind)
                        break;
                    end

                    backtrackCount = backtrackCount + 1;
                    mu0(:, ind) = 0.2 * mu0(:, ind); % shinkage
                    A_pr(:, ind) = A_pr(:, ind) + mu0(:, ind) .* searchDir(:, ind);
                    curr_y = Xfun(A_pr);
                    f1 = 0.5/p * vecnorm(abs(curr_y(:, ind))-ymag(:, ind), 2, 1).^2;
                end
            end
            
            if strcmp(opts.updateMu, 'True') 
                Dx = mu0 .* searchDir;
            end
            
            if strcmp(searchMethod, 'LBFGS')
                Dx = mu1 .* searchDir;
            end

            % save records
            if saveErr  % least-square measurement error
                errCurve(compt) = sum(vecnorm(ymag - abs(curr_y), 2, 1)); 
            end   
            if saveCorr % correlation coefficient to the true measurement
                corrCurve(compt) = sum(vecCorr(ymag, abs(curr_y), 1)); 
            end
            ts(compt) = toc(t0);
            compt = compt+1;
        end 

    elseif isfield(opts, 'times')
        errCurve = [];
        corrCurve = [];
        ts = [0];

        t0 = tic;
        compt = 1; compt0 = 1;
        
        while ts(end) <= opts.times        
            curr_y = Xfun(A_pr);
            weight = abs(curr_y) ./ (abs(curr_y)+alp*ymag); % gradient weight
            grad = Xinvfun((curr_y - ymag.*curr_y./abs(curr_y)).*weight)/p;
            
            % update via search direction & stepsize
            switch searchMethod
                case 'SD'
                    searchDir = -grad;
                case 'NCG'
                    if (compt==1 || compt-compt0 == CGresetPeriod) 
                        searchDir0 = zeros(size(grad), 'single');
                        beta = zeros(1, m, 'single');
                        compt0 = compt;
                    else       
                        if strcmpi(opts.beta, 'HS')
                            Dg = grad - grad0;
                            beta = -real(sum(conj(grad).*Dg, 1)) ./ real(sum(conj(searchDir0).*Dg, 1));
                        elseif strcmpi(opts.beta, 'FR')
                            beta = sum(grad.^2, 1) ./ sum(grad0.^2, 1);
                        end
                    end
                    searchDir = -grad + beta .*searchDir0;

                case 'LBFGS'
                    searchDir = -grad;
                    iterSub = min(compt-compt0, opts.stores);

                    if (compt==1 || compt-compt0 == BFGSresetPeriod) 
                        % Perform LBFGS initialization
                        yVals = zeros(n, m, opts.stores, 'single');
                        sVals = zeros(n, m, opts.stores, 'single');
                        rhoVals = zeros(1, m, opts.stores, 'single');
                        compt0 = compt;
                    else
                        % Update LBFGS stored vectors
                        Dg = grad - grad0; 
                        sVals = cat(3, Dx, sVals(:, :, 1:opts.stores-1));
                        yVals = cat(3, Dg, yVals(:, :, 1:opts.stores-1));
                        rhoVals = cat(3, 1./real(sum(conj(Dg).*Dx, 1)), rhoVals(:, :, 1:opts.stores-1));
                    end  
                    
                    if iterSub > 0  % two loops to determine the search direction
                        alphas = zeros(1, m, iterSub, 'single');
                
                        % First loop
                        for j = 1 : iterSub
                            alphas(:,:,j) = rhoVals(:,:,j) .* real(sum(conj(sVals(:,:,j)) .* searchDir, 1));
                            searchDir = searchDir - alphas(:,:,j) .* yVals(:,:,j);
                        end
                        
                        % Scaling of search direction
                        gamma = real(sum(conj(Dg) .* Dx, 1)) ./ sum(conj(Dg) .* Dg, 1);
                        searchDir = gamma .* searchDir;
                        
                        % Second loop
                        for j = iterSub : -1 : 1
                            beta = rhoVals(:,:,j) .* real(sum(conj(yVals(:,:,j)) .* searchDir, 1));
                            searchDir = searchDir + (alphas(:,:,j) - beta) .* sVals(:,:,j);
                        end
                        
                        searchDir = 1./gamma .* searchDir;
                        searchDir = vecnorm(grad, 2, 1) ./ vecnorm(searchDir, 2, 1) .* searchDir;
                    end
            end
            if strcmp(opts.updateMu, 'True') || strcmp(searchMethod, 'NCG') || strcmp(searchMethod, 'LBFGS')
                grad0 = grad;
            end
            if strcmp(opts.updateMu, 'True') || strcmp(searchMethod, 'NCG')
                searchDir0 = searchDir;
            end
 
            % update the stepsize with Barzilai-Borwein adaptive method
            if strcmp(opts.updateMu, 'True')
                if mod(compt-1, 3)==0                           % reinitialization of stepsize
                    mu1 = muRAF;
                else                                           % update the stepsize
                    Ds = searchDir - searchDir0;
                    dotprod = real(sum(conj(Dx).*Ds, 1));
                    muS = vecnorm(Dx,2,1).^2 ./ dotprod;       % First BB stepsize rule
                    muM = dotprod ./ vecnorm(Ds,2,1).^2;       % Alternate BB stepsize rule
                    muM = max(muM, 0);
                    ind = find(2*muM <= muS);
                    mu1 = muM; mu1(ind) = muS(ind)-muM(ind)./2; 
                    ind = union(find(mu1<=0), find(isinf(mu1))); 
                    ind = union(ind, find(isnan(mu1)));
                    if numel(mu0) > 1
                        mu1(ind) = 1.5*mu0(ind);              % grows first, then shinks with backtracking
                    else
                        mu1(ind) = 1.5*mu0;
                    end
                end
                f0 = 0.5/p * vecnorm(abs(curr_y)-ymag, 2, 1).^2;
            else
                mu1 = muRAF;
            end
            

            % update the TM
            A_pr = A_pr + mu1 .* searchDir;
            curr_y = Xfun(A_pr);
            
            if strcmp(opts.updateMu, 'True')
                mu0 = mu1; 
            end
 
            % Backtracking to find the stepsize with Armijo-Goldstein condition
            if strcmp(opts.updateMu, 'True') && mod(compt-1, 3)~=0
                ind = 1:m;
                f1 = 0.5/p * vecnorm(abs(curr_y(:, ind))-ymag(:, ind), 2, 1).^2;
                backtrackCount = 0;
                while backtrackCount <= 20 
                    tmp = f0(:, ind) + 0.1 * mu0(:, ind) .* real(sum(conj(searchDir0(:, ind)).*grad0(:, ind), 1));  
                    ind = ind(f1>tmp);
                    if isempty(ind)
                        break;
                    end

                    backtrackCount = backtrackCount + 1;
                    mu0(:, ind) = 0.2 * mu0(:, ind); % shinkage
                    A_pr(:, ind) = A_pr(:, ind) + mu0(:, ind) .* searchDir(:, ind);
                    curr_y = Xfun(A_pr);
                    f1 = 0.5/p * vecnorm(abs(curr_y(:, ind))-ymag(:, ind), 2, 1).^2;
                end
            end
            
            if strcmp(opts.updateMu, 'True') 
                Dx = mu0 .* searchDir;
            end
            
            if strcmp(searchMethod, 'LBFGS')
                Dx = mu1 .* searchDir;
            end

            % save records
            if saveErr  % least-square measurement error
                errCurve = [errCurve, sum(vecnorm(ymag - abs(curr_y), 2, 1))]; 
            end   
            if saveCorr % correlation coefficient to the true measurement
                corrCurve = [corrCurve, sum(vecCorr(ymag, abs(curr_y), 1))]; 
            end
            ts = [ts, toc(t0)];
            compt = compt+1;
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