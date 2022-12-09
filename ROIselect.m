function [dataROI] = ROIselect(data,ROI,sc)
% data is 3D dataset, selecting ROI in in-plane dimension

dataROI = sc*mean(mean(data( ROI(2)-round((ROI(4)-1)/2):ROI(2)+round((ROI(4)-1)/2) , ROI(1)-round((ROI(3)-1)/2):ROI(1)+round((ROI(3)-1)/2) , :) , 1) , 2);
end