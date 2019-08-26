function [rho,theta,houghSpace] = houghTransform_fast(theImage,rhoSampleFrequency, thetaSampleFrequency)

    %Define the hough space
    [width,height] = size(theImage);

    rhoLimit = round(norm([width height])-1);
    rho = (-rhoLimit:rhoSampleFrequency:rhoLimit);          
    theta = (-90:thetaSampleFrequency:89);

    numThetas = numel(theta);
    houghSpace = zeros(numel(rho),numThetas);

    %Find the "edge" pixels
    [yIndicies,xIndicies] = find(theImage);


    for i = (1:numThetas)
        rhos = xIndicies*cosd(theta(i)) + yIndicies*sind(theta(i)) ;
        houghSpace(:,i) = hist(rhos,rho); 
    end
    

    
    figure, imshow(imadjust(mat2gray(houghSpace)),[],'XData',theta,'YData',rho, 'InitialMagnification', 'fit');
    xlabel('\theta (degreees)'), ylabel('\rho');
    axis on, axis normal, hold on;
    colormap(hot)

end