function [] = plot_spectra(csi, x_start, y_start, x_size, y_size, f_first, f_end, figno_1, figno_2)
% csi : input spectra
% x_start, y_start : upper left corner of ROI
% x_size, y_size : size of ROI to show
% f_first, f_end : indices of first and last freqs to show, optional
% figno_1, figno_2 : figure numbers, optional


if nargin < 6
    f_first = 1;
    f_end = size(csi,3);
end

if nargin < 8
    figno_1 = 1;
    figno_2 = 2;
end


csi_ = abs(csi(x_start:x_start+x_size-1, y_start:y_start+y_size-1, f_first:f_end));
y_max = max(csi_(:));


pixels = size(csi,2) * size(csi,1);
nx = sqrt(pixels);


nnx = x_size;
nny = y_size;

widthx = 1/nnx*.9;
widthy = 1/nny*.9;
marginx = 1/nnx*0.05;
marginy = 1/nny*0.05;
originx = marginx:1/nnx:1;
originy = 1-1/nny+marginy:-1/nny:0;


xorig = x_start;
yorig = y_start;


x = xorig:xorig+nnx-1;
y = yorig:yorig+nny-1;


figure(figno_1);
for i=1:nnx,
    for j=1:nny,
%        subplot('Position', [originx(i), originy(j), widthx, widthy])
        subplot('Position', [originy(nny-j+1), originx(nnx-i+1), widthy, widthx])
        
        if nargin < 6
            plot(f_first:f_end, abs((squeeze(csi(x(i),y(j),:)))),'k'); 
                    axis([f_first,f_end, 0, y_max/1]); axis off;
        else
            plot(f_first:f_end, abs((squeeze(csi(x(i),y(j), f_first:f_end)))),'k', 'LineWidth', 2); 
                    axis([f_first,f_end, 0, y_max/1]); axis off;
        end

    end
end


ref_image = sum(abs(csi(:,:,f_first:f_end)),3);
ref_image = ref_image / max(abs(ref_image(:)));

% dB scale
% ref_image = 20*log10(sum(abs(csi(:,:,f_first:f_end)),3));

figure(figno_2); imagesc(abs(ref_image)); colormap jet; colorbar; axis square; hold on;
rectangle('Position',[y(1),x(1),y(end)-y(1),x(end)-x(1)],'EdgeColor','y', 'LineStyle','--');
hold off;

