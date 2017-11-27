function plotNoAvgTimeseries(metab_map2,metabs,vxls);
%%
%%
sz = size(metab_map2);
%%

metabCell = cell(length(metabs),sz(4));
%%
for i = 1:sz(4)
    for j = 1:length(metabs)
        metabCell{j,i} = squeeze(metab_map2(:,:,metabs(j),i));
        metabCell{j,i} = metabCell{j,i}(:);
    end
end
%%
inds = zeros(length(metabCell{1,1}),sz(4));
for i = 1:sz(4)
    inds(:,i) = metabCell{1,i} ~= 0;
end
ind = sum(inds,2); ind = ind == sz(4);
%%
horz = zeros(sum(ind),length(metabs),sz(4));
for i = 1:sz(4)
    for j = 1:length(metabs)
        horz(:,j,i) = metabCell{j,i}(ind);
    end
end
%%

lim = max(horz(:))*.9;
%%
% figure
% imagesc(cm);
% caxis([0 1]); 
% colormap jet; 
%%
if length(metabs) > 1
    figure('position',[0 1000 300*sz(4) 500])
    cum = 1;
    title('By run')
    for i = 1:sz(4);
            subplot(2,sz(4),cum)
            imagesc(metab_map2(:,:,metabs(1),i));
            axis square
            cum = cum + 1;

    end
    for i = 1:sz(4);
            subplot(2,sz(4),cum)
            imagesc(metab_map2(:,:,metabs(2),i));
            axis square
            cum = cum + 1;

    end
    tightfig;
    %%
    figure('position',[0 0 300*sz(4) 500])
    cum = 1;
    title('By run')
    for i = 1:sz(4);
            subplot(1,sz(4),cum)

            cm = corrcoef(horz(:,:,i));
            s = scatter(horz(:,1,i),horz(:,2,i),'filled');
            xl = xlim; yl = ylim;
            text(xl(1)+0.05,yl(2),...
                strcat('r = ',round(num2str(cm(1,2),2))),'fontsize',12);
            axis square
            alpha(0.2)
            cum = cum + 1;
    xlabel(['Conc metab: ' num2str(metabs(1))],'fontsize',16)
    ylabel(['Conc metab: ' num2str(metabs(2))],'fontsize',16)

    end
    tightfig;
    %%
    temp = horz(:,:,1);
    for i = 2:sz(4)
        temp = vertcat(temp,horz(:,:,i));
    end
    figure
    scatter(temp(:,1),temp(:,2),'filled')
    alpha(0.2)
    [r,p] = corr(temp);
    text(min(temp(:,1)),min(temp(:,2)),...
        strcat('r = ',round(num2str(r(1,2),2)),'; p = ', round(num2str(p(1,2),2))),...
        'fontsize',16);
    xlabel(['Conc metab: ' num2str(metabs(1))],'fontsize',16)
    ylabel(['Conc metab: ' num2str(metabs(2))],'fontsize',16)
    %%
    if size(vxls,1) > 0
        ss = get(0,'screensize');
        figure('position',[ss(3)-300 ss(4) 400 300*size(vxls,1)])
        for i = 1:size(vxls,1);

            subplot(size(vxls,1),1,i)
            hold on
            plot(squeeze(metab_map2(vxls(i,2),vxls(i,1),metabs(1),:)),'o-',...
                'linewidth',2)
            yyaxis right
            plot(squeeze(metab_map2(vxls(i,2),vxls(i,1),metabs(2),:)),'o-',...
                'linewidth',2)
            set(gca,'xtick',1:sz(4))
        end
    end
else
    figure('position',[0 0 300*sz(3) 400])
    for i = 1:sz(4)
        subplot(1,sz(4),i)
        imagesc(metab_map2(:,:,metabs,i))
        axis square
        caxis([0 lim])
    end
    tightfig;
    for i = 1:size(vxls,1);
        figure
        subplot(size(vxls,1),1,i)
        
        hold on
        plot(squeeze(metab_map2(vxls(i,2),vxls(i,1),metabs(1),:)),'o-')
    end
    
end
tightfig;
