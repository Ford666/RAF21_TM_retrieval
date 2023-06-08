% Compare the performances of different TM retrieval algorithms under the condition of 
% (1) different input sizes N
% (2) different sampling rates P/N
% (3) different signal-to-noise ratio  

% 0: PBR distribution at default settings, 1: PBR progression with time, 
% 2: PBR progression with sampling rates, 3: PBR progression with input size, 4: SNR

% Author: Shengfu Cheng
% Updated Date: 09/03/2023

addpath('./utils/');
addpath('./algs/');
close all


%% Parameters
M = 140; 
algName = ["prVBEM", "GGS2-1", "RAF2-1"];
N = 32;                                                                    % default input size
gamma = 5;                                                                 % default sampling rate
alpNoise = 0;                                                              % default noise level
N_list = [12, 16, 20, 24, 28, 32, 36, 40];                                 % input size
iter_list = [80, 160, 240, 320, 400, 600];                                 % iteration number
gamma_list = [2, 3, 4, 5, 6, 7, 8];                                        % sampling rates
alpNoise_list = [3e-3, 1e-2, 3.5e-2, 1e-1, 2e-1, 3.5e-1, 6e-1, 1, 2];      % SNR 

expType = 0; 
testNum = 1;

P = N^2 * gamma; 

optsVBEM = struct;
optsGGS21 = struct;
optsRAF21 = struct; optsRAF21.mu = 3;


%indices of the foci and background area (for PBR calculation)
Nscan = 20; 
focNum = Nscan^2;
Xidx = linspace(1, M, Nscan+4); Xidx = round(Xidx(3:end-2)); Yidx = Xidx;
[Xloc, Yloc] = meshgrid(Xidx, Yidx);
focInd = sub2ind([M, M], Yloc, Xloc); 

A = zeros(M, M, 'single'); A(M/2+10:M-10, M/2+10:M-10) = 1; bkgInd = find(A);

switch expType
    case 0 
        %% PBR distribution at default settings
        fprintf('Exp%d| Focus PBR distribution when N=%d, running for 20 s\n', expType, N);
        optsVBEM.times = 20;
        optsGGS21.times = 20;
        optsRAF21.times = 20;

        optsGGS21.ratio = 2/3;

        simTM = generate_tm(M^2, N^2);% simulated TM
        Ibkg = mean(abs(simTM * exp(1i*zeros(N^2, 1))).^2);


        SLMmsk = exp(1i*2*pi*rand(N^2, P, 'single'));
        SpekAmp_0 = abs(simTM * SLMmsk); SpekAmp_0 = SpekAmp_0 ./ max(max(SpekAmp_0));
        SpekAmp_noise = SpekAmp_0 + alpNoise*SpekAmp_0.*(rand(M^2, P, 'single')-0.5); %add uniform noise   
        subSpekAmp = SpekAmp_noise(focInd, :);
        
        %prVBEM
        [TM_vbem_sub, ts_vbem, errCurve_vbem] = prVBEM_tm_cpu(subSpekAmp, SLMmsk, optsVBEM);
              
        %GGS2-1
        [TM_ggs21_sub, ts_ggs21, errCurve_ggs21] = GGS21_tm_cpu(subSpekAmp, SLMmsk, optsGGS21);
        
        %RAF2-1
        [TM_raf21_sub, ts_raf21, errCurve_raf21] = RAF21_tm_cpu(subSpekAmp.^2, SLMmsk, optsRAF21);
        
  
        % Results comparison
        figure('color', 'w', 'position', [150 150 620 480]); %error curve
        maxErr = max([max(errCurve_vbem), max(errCurve_ggs21), max(errCurve_raf21)]);
        semilogy(ts_vbem(2:end), errCurve_vbem./maxErr, 'k-', ts_ggs21(2:end), errCurve_ggs21./maxErr, 'b-', ts_raf21(2:end), errCurve_raf21./maxErr,'r-', 'linewidth', 2);
        legend('prVBEM', 'GGS2-1', 'RAF2-1','fontname', 'Times New Roman Bold','fontsize', 12);
        ylim([0 1]); grid on; set(gca,'FontSize',18, 'LineWidth', 1); 
        xlabel(sprintf('t/s (N=%d, \\gamma=%d)', N^2, gamma), 'fontname', 'Times New Roman Bold', 'fontsize', 21); ylabel('Normalized focus PBR', 'fontname', 'Times New Roman Bold', 'fontsize', 21);

  
        focPatn_ref = simTM(focInd(:), :)'; focPatn_ggs21 = TM_ggs21_sub'; focPatn_vbem = TM_vbem_sub'; focPatn_raf21 = TM_raf21_sub'; 
        Ifoc_ref = abs(simTM * exp(1i*angle(focPatn_ref))).^2; Ifoc_ref_max = reshape(max(Ifoc_ref, [], 2), M, M);
        Ifoc_ggs21 = abs(simTM * exp(1i*angle(focPatn_ggs21))).^2; Ifoc_ggs21_max = reshape(max(Ifoc_ggs21, [], 2), M, M);
        Ifoc_vbem = abs(simTM * exp(1i*angle(focPatn_vbem))).^2; Ifoc_vbem_max = reshape(max(Ifoc_vbem, [], 2), M, M);
        Ifoc_raf21 = abs(simTM * exp(1i*angle(focPatn_raf21))).^2; Ifoc_raf21_max = reshape(max(Ifoc_raf21, [], 2), M, M);

        
        pbr_ref = max(Ifoc_ref,[], 1)./Ibkg; pbr_vbem = max(Ifoc_vbem,[], 1)./Ibkg;   
        pbr_ggs21 = max(Ifoc_ggs21,[], 1)./Ibkg; pbr_raf21 = max(Ifoc_raf21,[], 1)./Ibkg; 
        
        figure('color', 'w', 'position', [150 150 1250 550]); %focus scanning results
        subplot(241), imagesc(Ifoc_ref_max); title(sprintf('SimTM PBR=%.1f', mean(pbr_ref)), 'fontname', 'Times New Roman Bold','fontsize', 14),
        subplot(242), imagesc(Ifoc_vbem_max); title(sprintf('prVBEM PBR=%.1f', mean(pbr_vbem)), 'fontname', 'Times New Roman Bold','fontsize', 14),
        subplot(243), imagesc(Ifoc_ggs21_max); title(sprintf('GGS2-1 PBR=%.1f', mean(pbr_ggs21)), 'fontname', 'Times New Roman Bold','fontsize', 14),
        subplot(244), imagesc(Ifoc_raf21_max); title(sprintf('RAF2-1 PBR=%.1f', mean(pbr_raf21)), 'fontname', 'Times New Roman Bold','fontsize', 14),
        colormap(gray)

        subplot(245), histogram(pbr_ref, linspace(0,850, 10)); xlabel('Focus PBR', 'fontname', 'Times New Roman Bold','fontsize', 14), ylabel('Count', 'fontname', 'Times New Roman Bold','fontsize', 14);
        subplot(246), histogram(pbr_vbem, linspace(0,850, 10)); xlabel('Focus PBR', 'fontname', 'Times New Roman Bold','fontsize', 14), ylabel('Count', 'fontname', 'Times New Roman Bold','fontsize', 14);
        subplot(247), histogram(pbr_ggs21, linspace(0,850, 10)); xlabel('Focus PBR', 'fontname', 'Times New Roman Bold','fontsize', 14), ylabel('Count', 'fontname', 'Times New Roman Bold','fontsize', 14);
        subplot(248), histogram(pbr_raf21, linspace(0,850, 10)); xlabel('Focus PBR', 'fontname', 'Times New Roman Bold','fontsize', 14), ylabel('Count', 'fontname', 'Times New Roman Bold','fontsize', 14);

    
    case 1 
        %% PBR progression with time  
        fprintf('Exp%d| PBR progression with time \n', expType);
        times_list = [1, 2, 4, 8, 16, 30, 40, 50, 60];  % running time list (s)
        
        pbr_ref = pi/4*(N^2-1)+1;
        pbr_vbem_list = zeros(length(times_list), testNum); pbr_ggs21_list = zeros(length(times_list), testNum); pbr_raf21_list = zeros(length(times_list), testNum);    
        iters_vbem_list = zeros(length(times_list), testNum); iters_ggs21_list = zeros(length(times_list), testNum); iters_raf21_list = zeros(length(times_list), testNum);    
         
        for testInd = 1:testNum
            simTM = generate_tm(M^2, N^2); % simulated TM
            Ibkg = mean(abs(simTM * exp(1i*zeros(N^2, 1))).^2);

            P = gamma * N^2;
            SLMmsk = exp(1i*2*pi*rand(N^2, P, 'single'));
            SpekAmp = abs(simTM * SLMmsk); SpekAmp = SpekAmp ./ max(max(SpekAmp));
            subSpekAmp = SpekAmp(focInd, :);

            itInd = 1;
            for times = times_list
                optsVBEM.times = times;
                optsGGS21.times = times;
                optsRAF21.times = times;
                
                %prVBEM
                [TM_vbem_sub, ts_vbem, errCurve_vbem] = prVBEM_tm_cpu(subSpekAmp, SLMmsk, optsVBEM);

                %GGS2-1
                [TM_ggs21_sub, ts_ggs21, errCurve_ggs21] = GGS21_tm_cpu(subSpekAmp, SLMmsk, optsGGS21);
                

                %RAF2-1
                [TM_raf21_sub, ts_raf21, errCurve_raf21] = RAF21_tm_cpu(subSpekAmp.^2, SLMmsk, optsRAF21);                
      
                % iterations taken
                iters_vbem_list(itInd, testInd) = numel(ts_vbem); iters_ggs21_list(itInd, testInd) = numel(ts_ggs21); iters_raf21_list(itInd, testInd) = numel(ts_raf21); 
                
                % achievable PBR
                focPatn_vbem = TM_vbem_sub'; focPatn_ggs21 = TM_ggs21_sub'; focPatn_raf21 = TM_raf21_sub';
                Ifoc_vbem = abs(simTM * exp(1i*angle(focPatn_vbem))).^2; pbr_vbem_list(itInd, testInd) = mean(max(Ifoc_vbem,[], 1)./Ibkg);
                Ifoc_ggs21 = abs(simTM * exp(1i*angle(focPatn_ggs21))).^2; pbr_ggs21_list(itInd, testInd) = mean(max(Ifoc_ggs21,[], 1)./Ibkg);
                Ifoc_raf21 = abs(simTM * exp(1i*angle(focPatn_raf21))).^2; pbr_raf21_list(itInd, testInd) = mean(max(Ifoc_raf21,[], 1)./Ibkg);  
                
                itInd = itInd+1;
                
            end
            
            fprintf('Test %d|%d done\n', testInd, testNum);

        end

        % results
        fig = figure('color', 'w', 'position', [150 150 1200 500]); %focus scanning results
        subplot(121),
        plot(times_list, repelem(pbr_ref./pbr_ref, 1, numel(times_list)), 'm--', 'linewidth', 2); hold on;
        errorbar(times_list, mean(pbr_vbem_list./pbr_ref,2), 0.5*std(pbr_vbem_list./pbr_ref,0,2), 'k-', 'linewidth', 2); 
        errorbar(times_list, mean(pbr_ggs21_list./pbr_ref,2), 0.5*std(pbr_ggs21_list./pbr_ref,0,2), 'b-', 'linewidth', 2); 
        errorbar(times_list, mean(pbr_raf21_list./pbr_ref,2), 0.5*std(pbr_raf21_list./pbr_ref,0,2), 'r-', 'linewidth', 2); hold off;
        ylim([0 1.05])
        legend('Theoretical', 'prVBEM', 'GGS2-1', 'RAF2-1', 'fontsize', 14, 'location', 'southeast');
        grid on; set(gca,'FontSize',18, 'LineWidth', 1); 
        xlabel(sprintf('t/s (N=%d, \\gamma=%d)', N^2, gamma), 'fontsize', 21); ylabel('Normalized focus PBR', 'fontsize', 21);


        subplot(122),
%         plot(times_list, mean(iters_vbem_list,2), 'k+-', 'linewidth', 2); hold on;
%         plot(times_list, mean(iters_ggs21_list,2), 'b+-', 'linewidth', 2); 
%         plot(times_list, mean(iters_raf21_list,2), 'r+-', 'linewidth', 2); hold off;
        errorbar(times_list, mean(iters_vbem_list,2), std(iters_vbem_list,0,2), 'k-', 'linewidth', 2); hold on;
        errorbar(times_list, mean(iters_ggs21_list,2), std(iters_ggs21_list,0,2), 'b+-', 'linewidth', 2); 
        errorbar(times_list, mean(iters_raf21_list,2), std(iters_raf21_list,0,2), 'r+-', 'linewidth', 2); hold off;


        legend('prVBEM', 'GGS2-1', 'RAF2-1', 'fontsize', 12, 'location', 'northwest');
        grid on; set(gca,'FontSize',18, 'LineWidth', 1); 
        xlabel(sprintf('t/s (N=%d, \\gamma=%d)', N^2, gamma), 'fontsize', 21); ylabel('Iterations', 'fontsize', 21); 


    case 2 
        %% PBR progression with sampling rate (gamma)  
        fprintf('Exp%d| PBR progression with sampling rate\n', expType);
        N = 24; 
        pbr_ref = pi/4*(N^2-1)+1;
        
        pbr_vbem_list = zeros(length(gamma_list), testNum); pbr_ggs21_list = zeros(length(gamma_list), testNum); pbr_raf21_list = zeros(length(gamma_list), testNum);    
       
        for testInd = 1:testNum
            simTM = generate_tm(M^2, N^2); % simulated TM
            Ibkg = mean(abs(simTM * exp(1i*zeros(N^2, 1))).^2);

            gaInd = 1;
            for gamma = gamma_list
                P = gamma * N^2;
                SLMmsk = exp(1i*2*pi*rand(N^2, P, 'single'));
                SpekAmp = abs(simTM * SLMmsk); SpekAmp = SpekAmp ./ max(max(SpekAmp));
                subSpekAmp = SpekAmp(focInd, :);

                optsVBEM.times = 20; %20s is used for comparison under N=576, gamma=4
                optsGGS21.times = 20;
                optsRAF21.times = 20;

                %prVBEM
                [TM_vbem_sub, ts_vbem] = prVBEM_tm_cpu(subSpekAmp, SLMmsk, optsVBEM);
                

                %GGS2-1
                [TM_ggs21_sub, ts_ggs21] = GGS21_tm_cpu(subSpekAmp, SLMmsk, optsGGS21);
                

                %RAF2-1
                [TM_raf21_sub, ts_raf21, errCurve_raf21] = RAF21_tm_cpu(subSpekAmp.^2, SLMmsk, optsRAF21);
                if any(isnan(errCurve_raf21))
                    [TM_raf21_sub, ts_raf21] = RAF_tm_cpu(subSpekAmp.^2, SLMmsk, optsRAF21);
                end
                
                % achievable PBR
                focPatn_vbem = TM_vbem_sub'; focPatn_ggs21 = TM_ggs21_sub'; focPatn_raf21 = TM_raf21_sub';
                Ifoc_vbem = abs(simTM * exp(1i*angle(focPatn_vbem))).^2; pbr_vbem_list(gaInd, testInd) = mean(max(Ifoc_vbem,[], 1)./Ibkg);
                Ifoc_ggs21 = abs(simTM * exp(1i*angle(focPatn_ggs21))).^2; pbr_ggs21_list(gaInd, testInd) = mean(max(Ifoc_ggs21,[], 1)./Ibkg);
                Ifoc_raf21 = abs(simTM * exp(1i*angle(focPatn_raf21))).^2; pbr_raf21_list(gaInd, testInd) = mean(max(Ifoc_raf21,[], 1)./Ibkg);  
                
                gaInd = gaInd+1;
                
            end
            
            fprintf('Test %d|%d done\n', testInd, testNum);

        end
      
        % achievable PBR comparison
        fig = figure('color', 'w', 'position', [150 150 600 500]); %focus scanning results
        plot(gamma_list, repelem(pbr_ref./pbr_ref, 1, numel(gamma_list)), 'm--', 'linewidth', 2); hold on;
        errorbar(gamma_list, mean(pbr_vbem_list./pbr_ref,2), std(pbr_vbem_list./pbr_ref,0,2), 'k-', 'linewidth', 2); 
        errorbar(gamma_list, mean(pbr_ggs21_list/pbr_ref,2), std(pbr_ggs21_list./pbr_ref,0,2), 'b-', 'linewidth', 2); 
        errorbar(gamma_list, mean(pbr_raf21_list./pbr_ref,2), std(pbr_raf21_list./pbr_ref,0,2), 'r-', 'linewidth', 2); hold off;
        xlim([1.5 8]), ylim([0 1.05]);
        legend('Theoretical', 'prVBEM', 'GGS2-1', 'RAF2-1', 'fontsize', 14, 'location', 'southeast');
        grid on; set(gca,'FontSize',18, 'LineWidth', 1);
        xlabel(sprintf('\\gamma (N = %d, t=20 s)', N^2), 'fontsize', 21); ylabel('Normalized focus PBR', 'fontsize', 21);

    case 3 
        %% PBR progression with input size (N)
        fprintf('Exp%d| PBR progression with input size\n', expType);
        gamma = 4;
        pbr_ref = pi/4*(N_list.^2-1)+1; 
        pbr_vbem_list = zeros(length(N_list), testNum); pbr_ggs21_list = zeros(length(N_list), testNum); pbr_raf21_list = zeros(length(N_list), testNum); 
        t_list = [4,7,14,21,33,48,65,86];

        for testInd = 1:testNum
            simTM_all = generate_tm(M^2, max(N_list)^2); % simulated TM

            NInd = 1;
            for N = N_list
                simTM = simTM_all(:, 1:N^2);
                Ibkg = mean(abs(simTM * exp(1i*zeros(N^2, 1))).^2);

                P = gamma * N^2; 
                SLMmsk = exp(1i*2*pi*rand(N^2, P, 'single'));
                SpekAmp = abs(simTM * SLMmsk); SpekAmp = SpekAmp ./ max(max(SpekAmp));
                subSpekAmp = SpekAmp(focInd, :);
                
                % increase running time while keeping the iterations of
                % prVBEM fixed (e.g., 90)
                optsVBEM.times = t_list(NInd); 
                optsGGS21.times = t_list(NInd); 
                optsRAF21.times = t_list(NInd); 

                %prVBEM
                [TM_vbem_sub, ts_vbem] = prVBEM_tm_cpu(subSpekAmp, SLMmsk, optsVBEM);
                

                %GGS2-1
                [TM_ggs21_sub, ts_ggs21] = GGS21_tm_cpu(subSpekAmp, SLMmsk, optsGGS21);
                

                %RAF2-1
                [TM_raf21_sub, ts_raf21] = RAF21_tm_cpu(subSpekAmp.^2, SLMmsk, optsRAF21);
                

                % achievable PBR
                focPatn_vbem = TM_vbem_sub'; focPatn_ggs21 = TM_ggs21_sub'; focPatn_raf21 = TM_raf21_sub';
                Ifoc_vbem = abs(simTM * exp(1i*angle(focPatn_vbem))).^2; pbr_vbem_list(NInd, testInd) = mean(max(Ifoc_vbem,[], 1)./Ibkg);
                Ifoc_ggs21 = abs(simTM * exp(1i*angle(focPatn_ggs21))).^2; pbr_ggs21_list(NInd, testInd) = mean(max(Ifoc_ggs21,[], 1)./Ibkg);
                Ifoc_raf21 = abs(simTM * exp(1i*angle(focPatn_raf21))).^2; pbr_raf21_list(NInd, testInd) = mean(max(Ifoc_raf21,[], 1)./Ibkg);
                
                NInd = NInd+1;
            end
            
            fprintf('Test %d|%d done\n', testInd, testNum);
        end
        
     
        % achievable PBR comparison
        fig = figure('color', 'w', 'position', [150 150 600 500]); %focus scanning results
        plot(N_list.^2, pbr_ref, 'm--', 'linewidth', 2); hold on;
        errorbar(N_list.^2, mean(pbr_vbem_list,2), std(pbr_vbem_list, 0, 2), 'k-', 'linewidth', 2); 
        errorbar(N_list.^2, mean(pbr_ggs21_list,2), std(pbr_ggs21_list, 0, 2), 'b-', 'linewidth', 2); 
        errorbar(N_list.^2, mean(pbr_raf21_list,2), std(pbr_raf21_list, 0, 2), 'r-', 'linewidth', 2); hold off;
        legend('Theoretical', 'prVBEM', 'GGS2-1', 'RAF2-1', 'fontsize', 12, 'location', 'northwest');
        xlim([100 N_list(end)^2]); 
        grid on; set(gca,'FontSize',18, 'LineWidth', 1); 
        xticks(N_list.^2); xlabel(sprintf('N (\\gamma = %d)', gamma),'fontsize', 21); ylabel('Normalized focus PBR', 'fontsize', 21);
       

    case 4 
        %% Focusing PBR and uniformity with SNR (alp)
        fprintf('PBR progression with SNR\n');
        N = 24; gamma = 4;

        pbr_ref = pi/4*(N^2-1)+1;
        pbr_vbem_list = zeros(length(alpNoise_list), testNum); pbr_ggs21_list = zeros(length(alpNoise_list), testNum); 
        pbr_raf_list = zeros(length(alpNoise_list), testNum); pbr_raf21_list = zeros(length(alpNoise_list), testNum); 
        unif_vbem_list = zeros(length(alpNoise_list), testNum); unif_ggs21_list = zeros(length(alpNoise_list), testNum); 
        unif_raf_list = zeros(length(alpNoise_list), testNum);  unif_raf21_list = zeros(length(alpNoise_list), testNum); 
        
        SNR_list = zeros(length(alpNoise_list), testNum);
        
        %indices of the background region (for PBR calculation)
        for testInd = 1:testNum
            simTM = generate_tm(M^2, N^2);% simulated TM 
            Ibkg = mean(abs(simTM * exp(1i*zeros(N^2, 1))).^2);

            P = gamma*N^2;
            SLMmsk = exp(1i*2*pi*rand(N^2, P, 'single'));
            SpekAmp_0 = abs(simTM * SLMmsk); 
            subSpekAmp_0 = SpekAmp_0(focInd, :);
            subSpekAmp_0 = subSpekAmp_0./max(subSpekAmp_0, [], 'all');
                       
            alpNInd = 1;
            for alpNoise = alpNoise_list(1:end)  
                subSpekAmp_noise = subSpekAmp_0 + alpNoise*subSpekAmp_0.*(rand(focNum, P, 'single')-0.5); %add uniform noise
                SNR_list(alpNInd, testInd) = mean(subSpekAmp_0, 'all')/std(subSpekAmp_0-subSpekAmp_noise, 1, 'all');
                subSpekAmp_noise = subSpekAmp_noise./max(subSpekAmp_noise, [], 'all');

                optsVBEM.times = 20; %20s is used for comparison under N=576, gamma=4
                optsGGS21.times = 20;
                optsRAF21.times = 20;

                %prVBEM
                TM_vbem_sub = prVBEM_tm_cpu(subSpekAmp_noise, SLMmsk, optsVBEM);
    
                %GGS2-1e
                TM_ggs21_sub = GGS21_tm_cpu(subSpekAmp_noise, SLMmsk, optsGGS21);

                %RAF
                TM_raf_sub = RAF_tm_cpu(subSpekAmp_noise.^2, SLMmsk, optsRAF21);

                %RAF21
                [TM_raf21_sub, ~, errCurve_raf21] = RAF21_tm_cpu(subSpekAmp_noise.^2, SLMmsk, optsRAF21);

                % achievable PBR
                focPatn_vbem = TM_vbem_sub'; focPatn_ggs21 = TM_ggs21_sub'; focPatn_raf = TM_raf_sub'; focPatn_raf21 = TM_raf21_sub';     
                Ifoc_vbem = abs(simTM * exp(1i*angle(focPatn_vbem))).^2; pbr_vbem_list(alpNInd, testInd) = mean(max(Ifoc_vbem,[], 1)./Ibkg);
                Ifoc_ggs21 = abs(simTM * exp(1i*angle(focPatn_ggs21))).^2; pbr_ggs21_list(alpNInd, testInd) = mean(max(Ifoc_ggs21,[], 1)./Ibkg);
                Ifoc_raf = abs(simTM * exp(1i*angle(focPatn_raf))).^2; pbr_raf_list(alpNInd, testInd) = mean(max(Ifoc_raf, [], 1)./Ibkg);
                Ifoc_raf21 = abs(simTM * exp(1i*angle(focPatn_raf21))).^2; pbr_raf21_list(alpNInd, testInd) = mean(max(Ifoc_raf21, [], 1)./Ibkg);
                Ifoc_vbem_max = max(Ifoc_vbem, [], 2); Ifoc_ggs21_max = max(Ifoc_ggs21, [], 2); Ifoc_raf_max = max(Ifoc_raf, [], 2); Ifoc_raf21_max = max(Ifoc_raf21, [], 2);
                unif_vbem_list(alpNInd, testInd) = 1-std(Ifoc_vbem_max(focInd), 0, 'all')/mean(Ifoc_vbem_max(focInd), 'all'); unif_ggs21_list(alpNInd, testInd) = 1-std(Ifoc_ggs21_max(focInd), 0, 'all')/mean(Ifoc_ggs21_max(focInd), 'all'); 
                unif_raf_list(alpNInd, testInd) = 1-std(Ifoc_raf_max(focInd), 0, 'all')/mean(Ifoc_raf_max(focInd), 'all'); unif_raf21_list(alpNInd, testInd) = 1-std(Ifoc_raf21_max(focInd), 0, 'all')/mean(Ifoc_raf21_max(focInd), 'all');

                alpNInd = alpNInd+1;

            end
            
            fprintf('Test %d|%d done\n', testInd, testNum);
        end

        % achievable focus PBR & uniformity comparison
        fig = figure('color', 'w', 'position', [150 150 1100 450]); 
        subplot(121),
        semilogx(mean(SNR_list, 2), repelem(pbr_ref./pbr_ref, 1, length(alpNoise_list)), 'm--', 'linewidth', 2); hold on;
        errorbar(mean(SNR_list, 2), mean(pbr_vbem_list./pbr_ref,2), 0.5*std(pbr_vbem_list./pbr_ref, 0, 2), 'k-', 'linewidth', 2); 
        errorbar(mean(SNR_list, 2), mean(pbr_ggs21_list./pbr_ref,2), 0.5*std(pbr_ggs21_list./pbr_ref, 0, 2), 'b-', 'linewidth', 2); 
        errorbar(mean(SNR_list, 2), mean(pbr_raf_list./pbr_ref,2), 0.5*std(pbr_raf_list./pbr_ref, 0, 2), 'color', '#A2142F', 'linewidth', 2);
        errorbar(mean(SNR_list, 2), mean(pbr_raf21_list./pbr_ref,2), 0.5*std(pbr_raf21_list./pbr_ref, 0, 2), 'r-', 'linewidth', 2); hold off;
        
        legend('Theoretical', 'prVBEM', 'GGS2-1', 'RAF', 'RAF2-1', 'fontsize', 11, 'location', 'southeast');
        grid on; set(gca,'FontSize',18, 'LineWidth', 1); 
        xlim([0.8 1.1e3]); ylim([0 1.05]); xticks([1, 10, 10^2, 10^3])

        xlabel(sprintf('SNR (N = %d, \\gamma = %d, t=20s)', N^2, gamma), 'fontsize', 21); ylabel('Normalized focus PBR', 'fontsize', 21);

        ax2 = subplot(122);
        unif_vbem_list2 = unif_vbem_list; unif_vbem_list2(end,:) = 0.5*unif_vbem_list2(end,:);
        errorbar(mean(SNR_list, 2), mean(unif_vbem_list2,2), std(unif_vbem_list2, 0, 2), 'k-', 'linewidth', 2); hold on
        errorbar(mean(SNR_list, 2), mean(unif_ggs21_list,2), std(unif_ggs21_list, 0, 2), 'b-', 'linewidth', 2); 
        errorbar(mean(SNR_list, 2), mean(unif_raf_list,2), std(unif_raf_list, 0, 2), 'color', '#A2142F', 'linewidth', 2);
        errorbar(mean(SNR_list, 2), mean(unif_raf21_list,2), std(unif_raf21_list, 0, 2), 'r-', 'linewidth', 2); hold off;
        ax2.XScale = 'log'; 
        %legend('prVBEM', 'GGS2-1', 'RAF', 'RAF2-1', 'fontsize', 12, 'location', 'southeast');
        grid on; set(gca,'FontSize',18, 'LineWidth', 1); 
        xlim([0.8 1.1e3]); xticks([1, 10, 10^2, 10^3])
        xlabel(sprintf('SNR (N = %d, \\gamma = %d, t=20s)', N^2, gamma), 'fontsize', 21); ylabel('Focus uniformity', 'fontsize', 21);

end

