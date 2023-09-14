function [A_pr, ts, errCurve] = prVBEM_tm(y, x, iters, useGPU)
    % Phase recovery algorithm for TM retrieval
    % "Phase recovery from a Bayesian point of view: the variational approach" by A. Drémeau and F. Krzakala
    % The code is adopted from https://github.com/sudarshannagesh90/phase_retrieval/blob/master/Algorithms/prVBEM_withnoise
    
    saveErr = (nargout > 2);

    [m, p] = size(y);
    [n, ~] = size(x);

    y = y';
    H = x * x';
    
    if useGPU 
        x = gpuArray(x');
    else 
        x = x';
    end
    xH = Tikinv(x);

    % Initialization
    %*******************
    % default parameters
    tmp_A  = gather(xH * y);
    var_a  = max(abs(tmp_A), [], 1).^2; 
    phase = exp(1i*2*pi*rand(p,m,'single'));
    A_pr = xH * (y.*phase);
    var_n      = var(y, 0, 1)*1e-6;
    flag_est_n = 'on';
    pas_est    = 0.1;
    flag_cv    = 'KL';

    var_A   = zeros(n, m, 'single');

    ybar = y;
    z = x' * ybar;
        
    var_A = var_a .* var_n ./ (var_n + repelem(diag(H), 1, m) * diag(var_a)); 
    ym = gather(x * A_pr);
    
    errCurve = zeros(iters, 1, 'single');
    % Iterative process
    %*******************
    % Convergence criteria
    ts = zeros(iters, 1);
    t0 = tic;
    compt = 1;
    while compt <= iters

        % Estimation of var_n  
        %*********************
        if strcmp(flag_est_n,'on') || strcmp(flag_est_n,'on_off') 
            var_n = 1/p * (sum(abs(ym).^2, 1) + diag(H)'*var_A + sum(abs(y).^2, 1) - 2*real(sum(conj(ybar).*ym, 1)));
            if ~isreal(var_n)
                compt
                error('var_n is not real')
            end
        end

        % Update of the q(A_i)=G(A0,var_A)
        %************************************
        for k = 1:n     
            val_tmp = z(k, :)-H(k,:)*A_pr+H(k,k)*A_pr(k, :);
            A_pr(k, :) = (var_a./(var_n + var_a*H(k,k))).*val_tmp;
            var_A(k, :) = gather(var_n .* var_a./(var_n+var_a*H(k,k)));   
        end 

        t = gather(conj(y).*(x*A_pr));
        phi = t./abs(t);
        I0 = besseli(zeros(p,m,'single'), (2./var_n).*abs(t));
        I1 = besseli(ones(p,m,'single'), (2./var_n).*abs(t));
        fac_bessel = I1./I0;
        if ~isempty(find(isnan(fac_bessel)==1,1))
            fac_bessel(isnan(I1./I0)==1)=1;
        end
        ybar = y.*phi.*fac_bessel;
        z = x' * ybar;
        ym = gather(x * A_pr);    
        if saveErr; errCurve(compt) = gather(mean(vecnorm(y - abs(ym), 2, 1))); end
            
        ts(compt) = toc(t0);
        compt = compt+1;
    end
    A_pr = gather(A_pr');

    % normalization to std=1
    infInd = find(abs(A_pr)>1e10);
    A_pr(infInd) = exp(1i*angle(A_pr(infInd)));
    A_pr = 1/std(A_pr, 1, 'all') * A_pr;
    
end
