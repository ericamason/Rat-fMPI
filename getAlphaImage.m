function [Alp] = getAlphaImage(Image,Threshold,MaxOpacity,UpperBound)   
% Threshold: higher values show less of the overlay, lower values show more overlay
% MaxOpacit: This value is the maximum level of opacity in the overlay. a value of 1 indicates that the overlay can be fully opaque, .5 means it can be at most 50% opaque etc.
% UpperBound:

    Image = cast(Image,'double');
    Alp = ones(size(Image)); % Alp all ones;
    Alp((Image/max(Image,[],'all'))<Threshold) = 0; % any voxels less than threshold (% of max signal), set Alp = 0; 

    Image(Alp==0) = 0; % if below threshold, set to 0 in Image
    Image = Image-min(Image(Image>0),[],'all');
    Image(Image > UpperBound*max(Image,[],'all')) = UpperBound*max(Image,[],'all'); % EEM addition, cap all values above upperbound to that value (limits caxis)
    Image = Image/max(Image,[],'all')*255;

    Image(Image<0)=0; %making sure there are no negative values
%     Image2D = Image2D.^(1/2); %Compressing the dynamic range of the image...
    %Compressing the dynamic range reduces the speed that the transparency
    %falls off

    Alp = Image/max(Image,[],'all');
    Alp(Alp>MaxOpacity)=MaxOpacity;
    
end
