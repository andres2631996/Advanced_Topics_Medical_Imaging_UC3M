function [rho,theta,houghSpace] = HoughTransform_basic(theImage,rhoSampleFrequency, thetaSampleFrequency)

    %Define the hough space
    [width,height] = size(theImage);

    rhoLimit = round(norm([width height])-1);
    rho = (-rhoLimit:rhoSampleFrequency:rhoLimit);          
    theta = (-90:thetaSampleFrequency:89);

    numThetas = numel(theta);
    houghSpace = zeros(numel(rho),numThetas);

    %Find the "edge" pixels

    [yIndicies,xIndicies] = find(theImage);

    % Now you need a loop that iterates through all thetas
    % For every theta write the code to calculate the rho value for every
    % (x,y) pair resulting from the previous operation (edge pixels)
    % Update the houghSpace for those values of theta and rho
    for i=1:length(theta)
        for j=1:length(xIndicies)
            rho_initial=xIndicies(j)*cosd(theta(i))+yIndicies(j)*sind(theta(i));
            [aux,rho_pos]=min(abs(rho-rho_initial));
            houghSpace(rho_pos,i)=houghSpace(rho_pos,i)+1;
        end
    end
    

    % HoughSpace is ploted. It should look similar to the one you obtain
    % with Hough function in Matlab
    
    figure, imshow(imadjust(mat2gray(houghSpace)),[],'XData',theta,'YData',rho, 'InitialMagnification', 'fit');
    xlabel('\theta (degreees)'), ylabel('\rho');
    axis on, axis normal, hold on;
    colormap(hot)


end