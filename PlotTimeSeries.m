function [] = PlotTimeSeries(DataIn_conv, DriftTerms_All_RS, Yhat_RS, ConstReg_RS, CNRMap, ROI, TimeSec, CapniaTrigSeries_delay, conv_kernel_size, dx, dy, x, y, coreg_flag, coreg, bed_position, rows, a, flag_plotPSC, PSC, parench_mask)

fs_axislabels = 14;
fs_ticklabels = 12;
fs_scalebar = 14;
fs_titles = 14;

sc = 1e7;

%% Select ROI of time series
% DataIn_conv_ROI = sc*mean(mean(DataIn_conv( ROI(2)-round((ROI(4)-1)/2):ROI(2)+round((ROI(4)-1)/2) , ROI(1)-round((ROI(3)-1)/2):ROI(1)+round((ROI(3)-1)/2) , :) , 1) , 2);
% DriftTerms_All_RS_ROI = sc*mean(mean(DriftTerms_All_RS( ROI(2)-round((ROI(4)-1)/2):ROI(2)+round((ROI(4)-1)/2) , ROI(1)-round((ROI(3)-1)/2):ROI(1)+round((ROI(3)-1)/2) , :) , 1) , 2);
% Yhat_RS_ROI = sc*mean(mean(Yhat_RS( ROI(2)-round((ROI(4)-1)/2):ROI(2)+round((ROI(4)-1)/2) , ROI(1)-round((ROI(3)-1)/2):ROI(1)+round((ROI(3)-1)/2) , :) , 1) , 2);
% ConstReg_RS_ROI = sc*mean(mean(ConstReg_RS( ROI(2)-round((ROI(4)-1)/2):ROI(2)+round((ROI(4)-1)/2) , ROI(1)-round((ROI(3)-1)/2):ROI(1)+round((ROI(3)-1)/2)) , 1) , 2);

[DataIn_conv_ROI] = ROIselect(DataIn_conv,ROI,sc);
[DriftTerms_All_RS_ROI] = ROIselect(DriftTerms_All_RS,ROI,sc);
[Yhat_RS_ROI] = ROIselect(Yhat_RS,ROI,sc);
[ConstReg_RS_ROI] = ROIselect(ConstReg_RS,ROI,sc);

%% Plot time trace:
figure(201),
a1 = subplot(rows,4,[2:4]+4*(a-1));

% Plot yellow/orange activation boxes:
ylim([min(DataIn_conv_ROI-DriftTerms_All_RS_ROI) max(DataIn_conv_ROI-DriftTerms_All_RS_ROI)]); xlim([min(TimeSec) max(TimeSec)]),
ax = gca; ax.FontSize = fs_ticklabels;
[legVals] = plot_boxes_BlockActivation(CapniaTrigSeries_delay,TimeSec,0.9);

% Axis labels:
if a == rows
    xlabel('Time [sec]','fontweight','bold','fontsize',fs_axislabels,'Interpreter','latex'),
end
ylabel('ROI signal [A.U.]','fontweight','bold','fontsize',fs_axislabels-3,'Interpreter','latex');

% Plot data & regressor model fit:
hold on, plot(TimeSec,squeeze((DataIn_conv_ROI-DriftTerms_All_RS_ROI)) ,'LineWidth',1.5)
hold on, plot(TimeSec,squeeze((Yhat_RS_ROI-DriftTerms_All_RS_ROI)) ,'LineWidth',1.5)

% Add legend:
if a == 1
    f=(get(gca,'Children'));
    legend([f(1:2);f(end-1:end)],['Model fit','Measured data',legVals],'Location','northoutside','Interpreter','Latex')
end

% Add second Y-axis to show percent signal change:
a2 = axes('YAxisLocation','Right');
set(a2, 'color', 'none')
set(a2, 'XTick', [])
a2.FontSize = fs_ticklabels;
Perc = squeeze((DataIn_conv_ROI-DriftTerms_All_RS_ROI))*100/ConstReg_RS_ROI -100;
set(a2, 'YLim', [min(Perc) max(Perc)])
ylabel(a2,'$\Delta$ signal [\%]','fontweight','bold','fontsize',fs_axislabels,'Interpreter','latex');

% Link Y axes:
linkprop([a1,a2],'Position');

%% Plot CNR map
if coreg_flag
    [~,z_ind] = min(abs(coreg.z - bed_position));

    axH = subplot(rows,4,1+4*(a-1));
    Threshold = .1; MaxOpacity = 0.8; UpperBound = 0.3;
    %         Threshold = .1; MaxOpacity = 0.4; UpperBound = 0.3;
    [Alp] = getAlphaImage(CNRMap,Threshold,MaxOpacity,UpperBound); % can change Threshold, MaxOpacity, and/or UpperBound
    ax1 = axes; imagesc(ax1, x, y, coreg.MRI(:,:,z_ind)); set(gca,'YDir','normal'); axis image;
    ax2 = axes; imagesc(ax2, x, y, CNRMap,'AlphaData', Alp(:,:)); set(gca,'YDir','normal'); axis image; colorbar('location','westoutside', 'FontSize', fs_ticklabels);
    cl = caxis; caxis([0,cl(2)]);

    % Plot ROI:
    ROIx = [x(ROI(1)), x(ROI(1))+ROI(3)*dx];
    ROIy = [y(ROI(2)), y(ROI(2))+ROI(4)*dy];
    ROIx_ck = [min(ROIx) - abs(dx)*(conv_kernel_size-1)/2, max(ROIx) + abs(dx)*(conv_kernel_size-1)/2]; % with convolution kernel size added
    ROIy_ck = [min(ROIy) - abs(dy)*(conv_kernel_size-1)/2, max(ROIy) + abs(dy)*(conv_kernel_size-1)/2];
    hold on, plot(ax2, [min(ROIx_ck) max(ROIx_ck)],[min(ROIy_ck) min(ROIy_ck)],'c'); xlim([min(x) max(x)]);  ylim([min(y) max(y)]);
    hold on, plot(ax2, [min(ROIx_ck) max(ROIx_ck)],[max(ROIy_ck) max(ROIy_ck)],'c'); xlim([min(x) max(x)]);  ylim([min(y) max(y)]);
    hold on, plot(ax2, [min(ROIx_ck) min(ROIx_ck)],[min(ROIy_ck) max(ROIy_ck)],'c'); xlim([min(x) max(x)]);  ylim([min(y) max(y)]);
    hold on, plot(ax2, [max(ROIx_ck) max(ROIx_ck)],[min(ROIy_ck) max(ROIy_ck)],'c'); xlim([min(x) max(x)]);  ylim([min(y) max(y)]);

    % Plot spatial convolution kernel scale bar:
    hold on, plot(ax2, [x(4) x(4)+(dx*conv_kernel_size)],[y(end-4) y(end-4)],'Color','w','LineWidth',1)
    hold on, plot(ax2, [x(4) x(4)],[y(end-3) y(end-5)],'Color','w','LineWidth',1)
    hold on, plot(ax2, [x(4)+(dx*conv_kernel_size) x(4)+(dx*conv_kernel_size)],[y(end-3) y(end-5)],'Color','w','LineWidth',1)
    hold on, text(ax2, [x(8)+(dx*conv_kernel_size) x(8)+(dx*conv_kernel_size)] , [y(end-5) y(end-5)] , [num2str(dx*conv_kernel_size),' mm'] ,'Color','w','FontSize',fs_scalebar)

    linkprop([axH,ax1,ax2],'Position');  ax2.Visible = 'off'; axH.Visible = 'off'; ax1.XTick = []; ax1.YTick = []; %ax1.Visible = 'off';  ax2.YTick = [];
    colormap(ax1,'gray'); colormap(ax2,'redwhite');

    CNRMap_ROI = ROIselect(CNRMap,ROI,1);
    title(ax1,['CNR map',', max = ',num2str(roundn(CNRMap_ROI,-2))],'fontsize',fs_titles)
    % title(ax1,['CNR map',', max = ',num2str(roundn(max(CNRMap,[],'all'),-2))],'fontsize',fs_titles)

else
    ROIx = [ROI(1), ROI(1)+ROI(3)];
    ROIy = [ROI(2), ROI(2)+ROI(4)];

    subplot(rows,4,1+4*(a-1)),
    imagesc(CNRMap); axis image; colormap('redwhite'); colorbar('FontSize', fs_ticklabels)
    title(['CNR map',', ROI = ',num2str(roundn(mean(CNRMap(ROIy(1):ROIy(2), ROIx(1):ROIx(2)),'all'),-1))],'fontsize',fs_titles)
    hold on, rectangle('Position',ROI,'EdgeColor','k','LineWidth',1)
    caxis([-max(CNRMap(:)) max(CNRMap(:))]);
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);

end


%% plot PSC map:

 figure(202)
if flag_plotPSC
    % subplot(rows,1,a),
    % imagesc(PSC); axis image; colormap('parula'); colorbar('FontSize', fs_ticklabels)
    % set(gca,'YTickLabel',[]);
    % set(gca,'XTickLabel',[]);
    % clim([0 max(PSC,[],'all')]);
    % 
    % ROIx = [(ROI(1)), (ROI(1))+ROI(3)];
    % ROIy = [(ROI(2)), (ROI(2))+ROI(4)];
    % ROIx_ck = [min(ROIx) - (conv_kernel_size-1)/2, max(ROIx) + (conv_kernel_size-1)/2]; % with convolution kernel size added
    % ROIy_ck = [min(ROIy) - (conv_kernel_size-1)/2, max(ROIy) + (conv_kernel_size-1)/2];
    % hold on, plot([min(ROIx_ck) max(ROIx_ck)],[min(ROIy_ck) min(ROIy_ck)],'k'); 
    % hold on, plot([min(ROIx_ck) max(ROIx_ck)],[max(ROIy_ck) max(ROIy_ck)],'k'); 
    % hold on, plot([min(ROIx_ck) min(ROIx_ck)],[min(ROIy_ck) max(ROIy_ck)],'k'); 
    % hold on, plot([max(ROIx_ck) max(ROIx_ck)],[min(ROIy_ck) max(ROIy_ck)],'k'); 
    % PSC_ROI = ROIselect(PSC,ROI,1);
    % title(['dS/S [%]',', ROI = ',num2str(roundn(PSC_ROI,-1))],'fontsize',fs_titles)


    axH = subplot(rows,1,a); %  subplot(rows,4,1+4*(a-1));
    ax1 = axes; imagesc(ax1, x, y, PSC); set(gca,'YDir','normal'); axis image; colorbar('location','eastoutside', 'FontSize', fs_ticklabels);
    caxis([0 max(PSC,[],'all')])
    % caxis([0 12])
    outline_parench = edge(parench_mask); 
    ax2 = axes; imagesc(ax2, x, y, outline_parench,'AlphaData', outline_parench>0); set(gca,'YDir','normal'); axis image; 
    % cl = caxis; caxis([0,cl(2)]);

    % Plot ROI:
    % ROIx = [x(ROI(1)), x(ROI(1))+ROI(3)*dx];
    % ROIy = [y(ROI(2)), y(ROI(2))+ROI(4)*dy];
    % ROIx_ck = [min(ROIx) - abs(dx)*(conv_kernel_size-1)/2, max(ROIx) + abs(dx)*(conv_kernel_size-1)/2]; % with convolution kernel size added
    % ROIy_ck = [min(ROIy) - abs(dy)*(conv_kernel_size-1)/2, max(ROIy) + abs(dy)*(conv_kernel_size-1)/2];
    % hold on, plot(ax2, [min(ROIx_ck) max(ROIx_ck)],[min(ROIy_ck) min(ROIy_ck)],'c'); xlim([min(x) max(x)]);  ylim([min(y) max(y)]);
    % hold on, plot(ax2, [min(ROIx_ck) max(ROIx_ck)],[max(ROIy_ck) max(ROIy_ck)],'c'); xlim([min(x) max(x)]);  ylim([min(y) max(y)]);
    % hold on, plot(ax2, [min(ROIx_ck) min(ROIx_ck)],[min(ROIy_ck) max(ROIy_ck)],'c'); xlim([min(x) max(x)]);  ylim([min(y) max(y)]);
    % hold on, plot(ax2, [max(ROIx_ck) max(ROIx_ck)],[min(ROIy_ck) max(ROIy_ck)],'c'); xlim([min(x) max(x)]);  ylim([min(y) max(y)]);


    linkprop([axH,ax1,ax2],'Position');  ax2.Visible = 'off'; axH.Visible = 'off'; ax1.XTick = []; ax1.YTick = []; %ax1.Visible = 'off';  ax2.YTick = [];
    colormap(ax1,'parula'); colormap(ax2,'gray');

    PSCmask = PSC.*parench_mask; 
    PSCmask(PSCmask == 0) = nan; 
    mean_parench_PSC = mean(PSCmask, 'all','omitnan');
    % title(ax1,{['dS/S [%]'],['mean value in segmented region = ',num2str(roundn(mean_parench_PSC,-1)),'%']},'fontsize',12)
    title(ax1,['dS/S [%]'],'fontsize',fs_titles)
    % text(ax1, -max(x)*3/4, -max(y)/4, ['mean dS/S = ',num2str(roundn(mean_parench_PSC,-1)),'%'],'Color','w'); 

    %%%% 
    figure(203),
    axH = gca; 
    ax1 = axes; imagesc(ax1, x, y, coreg.MRI); set(gca,'YDir','normal'); axis image; colorbar('location','eastoutside', 'FontSize', fs_ticklabels);
    outline_parench = edge(parench_mask); 
    ax2 = axes; imagesc(ax2, x, y, outline_parench,'AlphaData', outline_parench>0); set(gca,'YDir','normal'); axis image; 

    linkprop([axH,ax1,ax2],'Position');  ax2.Visible = 'off'; axH.Visible = 'off'; ax1.XTick = []; ax1.YTick = []; %ax1.Visible = 'off';  ax2.YTick = [];
    colormap(ax1,'gray'); colormap(ax2,'gray');
end

end
