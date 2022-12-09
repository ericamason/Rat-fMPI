function [legend_order,pat] = plot_boxes_BlockActivation(CapniaTrigSeries,TimeSec, alpha)

if nargin < 3
    alpha = 0.9;
end

ChangeIndex = find(abs(diff(CapniaTrigSeries))==1);
ChangeIndex = [1 ChangeIndex length(CapniaTrigSeries)];

orange = [243 184 104]/255; %% HYPER
yellow = [254 237 121]/255; %% HYPO

for i = 1:(length(ChangeIndex)-1)
    CurrentYLimits = get(gca,'YLim');
    
    IsHyperCap = CapniaTrigSeries(ChangeIndex(i)+1)==1;
    if IsHyperCap
        pat(i) = patch([TimeSec(ChangeIndex(i)) TimeSec(ChangeIndex(i)) TimeSec(ChangeIndex(i+1)) TimeSec(ChangeIndex(i+1))], [CurrentYLimits(1) CurrentYLimits(2) CurrentYLimits(2) CurrentYLimits(1)],'k','FaceColor', orange , 'FaceAlpha', alpha);
        
    else
        pat(i) = patch([TimeSec(ChangeIndex(i)) TimeSec(ChangeIndex(i)) TimeSec(ChangeIndex(i+1)) TimeSec(ChangeIndex(i+1))], [CurrentYLimits(1) CurrentYLimits(2) CurrentYLimits(2) CurrentYLimits(1)],'k','FaceColor', yellow , 'FaceAlpha', alpha);
    end
    
    if i == 1
       if IsHyperCap
           legend_order = {'Hypercapnia','Hypocapnia'};
       else
           legend_order = {'Hypocapnia','Hypercapnia'};
       end
    end
end
% hLegend = findobj(gcf, 'Type', 'Legend');
% hLegend.String= {hLegend.String{1},hLegend.String{2}};
set(gca,'children',flipud(get(gca,'children')))
% xlim([0 TimeSec(end)])