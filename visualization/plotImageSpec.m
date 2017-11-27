function plotImageSpec(inputImgSpec, con, img2)

f1 = figure();
if con == 1
    inputImgSpec = conj(inputImgSpec);
    img2 = conj(img2);
end
figure(f1)
colormap jet
imagesc(sum(abs(inputImgSpec),3))
axis square
drawnow
[Vx,Vy] = ginput;
%%
Vx = floor(Vx); Vy = floor(Vy);
nnx = 3;
nny = 3;

widthx = 1/nnx*.9;
widthy = 1/nny*.9;
marginx = 1/nnx*0.05;
marginy = 1/nny*0.05;
originx = marginx:1/nnx:1;
originy = 1-1/nny+marginy:-1/nny:0;

%y_max = max(real(inputImgSpec(:)));
%y_min= min(real(inputImgSpec(:)));
    f_first = 512;
    f_end = 1024 %size(inputImgSpec,3);

for m = 1:length(Vx)
    
    figure('name',['Voxel coords : ' num2str(Vx(m)) ' ' num2str(Vy(m)) ]); 

    xorig = Vx(m);
    yorig = Vy(m);


    x = xorig-1:xorig+1;
    y = yorig-1:yorig+1;


    % x : vertical axis
    % y : horizontal axis
%% bug makes this not fill in the right order.
c = 0;
tot = nny*nnx;
    for j=1:nny,
        for i=1:nnx,
            c = c+1;
            subplot(3,3,c)%('Position', [originx(nnx-i+1), originy(nny-j+1), widthy, widthx])
                if isempty(img2)
                    hold on, plotSpecFinish(inputImgSpec,y(j),x(i),1)
                    xlim([f_first,f_end]); axis off;
                else
                   hold on, plotSpecFinish(inputImgSpec,y(j),x(i),1)
                   plotSpecFinish(img2,y(j),x(i),1),axis off
                end
                    
            
        end

    end
    tightfig;
    %%
    figure(f1)
    hold on
    rectangle('Position',[x(1),y(1),x(end)-x(1),y(end)-y(1)],'EdgeColor','k', 'LineStyle','--','Linewidth',2);

end


%%
