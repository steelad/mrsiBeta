function dispmatrix(A,slideryes)
cla
% % Define the matrix
% A = 999*rand(20,25);

% Determine the number of rows and columns
[m,n] = size(A);

% Draw the first column
axis([1 m 1 n])
drawnow
tmp = text(.5,.5,'t');
ext = get(tmp,'Extent');
dy = ext(4); % Height of a single character.
wch = ext(3); % Width of a single character.
fw = 8; % Column width. Each column can have up to
% 8-digits in it.
wc = 8*wch; % Width of 8 digit column
dwc = 2*wch; % Distance between columns
dx = wc+dwc; % Step used for columns
x = 1; % Location of first column
delete(tmp)
for i = 1:n % Column
y = 1;
for j = 1:m % Row
y = y + abs(dy); % Location of row
t(j,i) = text(x,y,sprintf('%3.4f',A(j,i)));
end
x = x+dx; % Location of next column
end

% Fix up the axes
axis([1-dwc/2 1+6*dx-dwc/2 1 n])
set(gca,'XTick',(1-dwc/2):dx:x)
set(gca,'XGrid','on','GridLineStyle','-','Box','on')
set(gca,'YTick',[],'XTickLabel',[],'Visible','on')
title('Columns 1 through 6')

if slideryes
% Add a horizontal slider
slide_step = floor(10/dx); % Number of columns in view
hs = uicontrol('Style','slider','Units','normal', ...
'Position',[.1 0 .8 .05],'min',1,'max',n, ...
'UserData',[dx,dwc],'Value',1, ...
'CallBack',['val = get(gco,''Value'');', ...
'dx = get(gco,''UserData'');', ...
'dwc = dx(2); dx = dx(1);', ...
'if (val-round(val))>eps,', ... % Test for + step
' val = ceil(val);', ...
'elseif (val-round(val)) < -eps,', ... % Test for - step
' val = floor(val);', ...
'end,', ...
'minx = 1+(floor(val)-1)*dx-dwc/2;', ...
'maxx = minx+6*dx;', ...
'set(gco,''Value'',val),', ... % Set the step
'axlim = axis;', ...
'axis([minx maxx axlim(3:4)]),', ...
'title([''Columns '',int2str(val),'' through '', int2str(val+6)])']);
end