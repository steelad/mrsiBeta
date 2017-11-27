function [] = overplot_spectra(csi1, csi2, x_start, y_start, x_size, y_size, f_first, f_end, figno_1, figno_2)
% csi : input spectra
% x_start, y_start : upper left corner of ROI
% x_size, y_size : size of ROI to show
% f_first, f_end : indices of first and last freqs to show, optional
% figno_1, figno_2 : figure numbers, optional

if nargin < 7
    f_first = 1;
    f_end = size(csi1,3);
end


csi_ = real(csi2(x_start:x_start+x_size-1, y_start:y_start+y_size-1, f_first:f_end));
y_max = max(csi_(:));
y_min= min(csi_(:));

pixels = size(csi1,2) * size(csi1,1);
nx = sqrt(pixels);
ny = nx;



if nargin < 9
   figno_1 = 1;
   figno_2 = 2;
end

h = figure(figno_1); close(h);
h = figure(figno_2); close(h);

   
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


% x : vertical axis
% y : horizontal axis


figure(figno_1); 
for j=1:nny,
    for i=1:nnx,
        subplot('Position', [originy(nny-j+1), originx(nnx-i+1), widthy, widthx])
                
        if nargin < 6
            hold on, plot(f_first:f_end, real((squeeze(csi1(x(i),y(j),:)))),'b'); 
            hold on, plot(f_first:f_end, real((squeeze(csi2(x(i),y(j),:)))),'k'); 
                axis([f_first,f_end, y_min, y_max / 1]); axis off;
        else
            hold on, plot(f_first:f_end, real((squeeze(csi1(x(i),y(j), f_first:f_end)))),'b', 'LineWidth', 2); 
            hold on, plot(f_first:f_end, real((squeeze(csi2(x(i),y(j), f_first:f_end)))),'k', 'LineWidth', 2); 
                axis([f_first,f_end, y_min, y_max / 1]); axis off;
        end

    end
end


ref_image = sum(abs(  cat(2, csi1, csi2)   ),3);

figure(figno_2);
imagesc( 20*log10(ref_image)); colormap jet; colorbar; axis square; hold on;
rectangle('Position',[y(1),x(1),y(end)-y(1),x(end)-x(1)],'EdgeColor','k', 'LineStyle','--');
rectangle('Position',[y(1)+ny,x(1),y(end)-y(1),x(end)-x(1)],'EdgeColor','k', 'LineStyle','--'); axis image
hold off;

