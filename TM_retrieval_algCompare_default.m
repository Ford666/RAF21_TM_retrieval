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


%% Focusing PBR distribution at default settings
fprintf('Focus PBR distribution when N=%d, running for 20 s\n', N);
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

print(gcf,'-dpng', '-r400', sprintf('Algorithm_comparison_errorCurve_default.png'));   



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

% print(gcf,'-dpng', '-r400', sprintf('Algorithm_comparison_image_histogram_default.png'));   
