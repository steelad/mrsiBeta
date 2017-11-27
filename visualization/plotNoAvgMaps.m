function plotNoAvgMaps(metab_map2,metab)
%%
sz = size(metab_map2);

%%
try
    metabCell = cell(1,sz(4));
catch
    sz(4) = 1;
end
%%
for i = 1:sz(4)
    metabCell{i} = squeeze(metab_map2(:,:,metab,i));
    metabCell{i} = metabCell{i}(:);
end
inds = zeros(length(metabCell{1}),sz(4));
for i = 1:sz(4)
    inds(:,i) = metabCell{i} ~= 0;
end
ind = sum(inds,2); ind = ind == sz(4);
%%
horz = zeros(sum(ind),sz(4));
for i = 1:sz(4)
    horz(:,i) = metabCell{i}(ind);
end
    lim = max(horz(:))*.9;
    cm = corrcoef(horz);
%%
% figure
% imagesc(cm);
% caxis([0 1]); 
% colormap jet; 
%%
figure('position',[0 0 1000 1000])
cum = 1;
for i = 1:sz(4);
    for j = 1:sz(4)
        if i ~= j
            subplot(sz(4),sz(4),cum)
            s = scatter(horz(:,i),horz(:,j),'filled');
            alpha(0.2)
            cum = cum + 1;
            text(min(horz(:,i)),max(horz(:,j)),...
                strcat('r = ',round(num2str(cm(i,j),2))));
            axis off
        else
            subplot(sz(4),sz(4),cum)
            histogram(horz(:,i),20,'normalization','probability')
            xlim([min(horz(:)) max(horz(:))])
            cum = cum + 1;
            
        end
    end
end
tightfig;
%%
%%
figure('position',[0 0 1200 300])
for i = 1:sz(4)
    subplot(1,sz(4),i)
    imagesc(metab_map2(:,:,metab,i))
    axis square
    caxis([0 lim])
end
tightfig;