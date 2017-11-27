
function [r,p,gm_wm_pv,gm_wm_t,meanCLRB] = corr_metab_gm(subjStruct,metab,metabName)

for i = 1:size(subjStruct.metab_map,4)
    ind = subjStruct.CLRB_map(:,:,metab,i) < 50 &...
        subjStruct.metab_map(:,:,metab,i) > 0 &...
        subjStruct.metab_map(:,:,metab,i) < 20;
    %%
    %%
    X = subjStruct.gm_prob(ind); 
    Y = subjStruct.metab_map(:,:,metab,i); Y = Y(ind);
    X = X(:); Y = Y(:);
    [r,p] = corr(X,Y);
    %%
    figure;
    scatter(X,Y,40,Y,'linewidth',3)
    hold on
    addLinFit(X,Y,[0.5 0.5 0.5],[0.2 max(Y)])
    title([subjStruct.subj ' ' metabName]);
    %ylim([0 15])
    colormap('parula')
    axis square
    formatFigAxis(gca,5,[],5,0,[min(X) max(Y)])
    %print(gcf,[subjStruct.subj '.' metabName '.scat.eps'],'-depsc')
    %close
    tempn = subjStruct.CLRB_map(:,:,metab);
    meanCLRB = mean(tempn(ind))
    %%
    %GM_Bin = X > 0.5;
    Y_above = Y(X > 0.9);
    Y_below = Y(X < 0.1);
     [~,gm_wm_pv,~,gm_wm_t] = ttest2(Y_above,Y_below);
    % %%
    % clr = lines(4);
    % scatJitt({Y_below,Y_above},clr(1:2,:),[],0)
    % title([subjStruct.subj ' ' metabName]);
    % formatFigAxis(gca,5,[1 2],{'WM','GM'},0,[])
    % text(0.55,min(Y_below),['t = ' sprintf('%.2f',gm_wm_t.tstat) ...
    %     '; p = ' sprintf('%.2f',gm_wm_pv)],'fontsize',12)
    % print(gcf,[subjStruct.subj '.' metabName '.split.eps'],'-depsc')
    % close
    %%
end