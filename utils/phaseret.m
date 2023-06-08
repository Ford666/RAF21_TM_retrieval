function [output] = phaseret(data,ref,mask)

    %% Spatial filtering
    scale = 1; 
    inters = single(data);
    [~,M,num] = size(inters); 
    
    t0=tic;
    Finters = fftshift(fftshift(fft2(inters, scale*M, scale*M), 1), 2); 
    fprintf("2D DFT takes %.2fs\n",toc(t0));

    Iref = im2double(ref);
    %select filter window manually
    Fh = abs(Finters(:, :, 1)).^2;
    if mask == 0
		figure('position', [200, 400, 550, 500]),
		imshow(Fh, [0,max(max(Fh))./1e3]), colormap(pink);

		%circle filter
		[x0, y0] = ginput(1);   %center
		[x1, y1] = ginput(1);   %circle edge
		[X, Y] = meshgrid(1:scale*M, 1:scale*M);
		r2 = (x0-x1)^2+(y0-y1)^2;
		Mask1 = (X-x0).^2 + (Y-y0).^2 < r2;

%    % rectangle filter mask
%     [x0, y0] = ginput(1);   %top-left
%     [x1, y1] = ginput(1);   %bottom-right
%     Mask1 = zeros(scale*M, scale*M, 'single');
%     Mask1(round(y0):round(y1), round(x0):round(x1)) = 1;

		save('./utils/phasefm.mat','Mask1');
    else

		load('./utils/phasefm.mat');
    end
    
    %spatial filtering for all interferences
    t0=tic;
    Finters_2 = Finters .* Mask1;

    clear Finters
    %zero-frequency centered and recover object
    Finters_3 = zeros(size(Finters_2), 'single');

    [row, col] = find(Mask1~=0); 

    minrow = round(min(row)); mincol = round(min(col)); 
    meanrow = round(mean(row)); meancol = round(mean(col)); 
    maxrow = round(max(row)); maxcol = round(max(col));

    Finters_3((minrow:maxrow)-meanrow+scale*M/2, (mincol:maxcol) - ...
            meancol+scale*M/2, :) = Finters_2(minrow:maxrow, mincol:maxcol, :);
    clear Finters_2
    fprintf("-1 order frequency shift takes %.2fs\n",toc(t0));

    t0=tic;
%     tempIFT = zeros(M*scale, M*scale, num, 'single');
%     N_worker = 2;
%     IdxBlock = num/N_worker;
%     for x = 1:N_worker
%         tempIFT(:, :, IdxBlock*(x-1)+1:IdxBlock*x) = ...
%                     ifft2(ifftshift(ifftshift(Finters_3(...
%                     :, :, IdxBlock*(x-1)+1:IdxBlock*x),1),2));
%     end
    tempIFT = ifft2(ifftshift(ifftshift(Finters_3,1),2));
    fprintf("2D IDFT takes %.2fs\n",toc(t0));

    t0=tic;
    output = tempIFT(1:M,1:M, :) ./sqrt(Iref);
%     output = tempIFT(1:M,1:M, :);
    clear Finters_3

    fprintf("Divided by reference intensity takes %.2fs\n",toc(t0));
end
