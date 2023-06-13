function himmelblauContours(edges)
% himmelblauContours   Generates a contour plot of the Himmelblau problem
%
% Syntax: himmelblauContours(edges)
%
%         NO OUTPUT

% Created: 06.06.2021
% Authors : Anton Ruby Larsen and Carl Frederik Grønvald
%           IMM, Technical University of Denmark

%%
    % How wide the plot should be
    if nargin<1
        edges = 5;
    end

    Himmelblau = @(x1,x2)(x1.^2+x2-11).^2 + (x1+x2.^2-7).^2;

    % Compute the function values
    x1 = linspace(-edges,edges, 400); % linspace(-5,5, 100);
    x2 = linspace(-edges,edges, 400); % linspace(-5,5, 100);
    [X1,X2]=meshgrid(x1,x2);
    F = Himmelblau(X1,X2);
    

    % Make contour plot
    if edges == 5
        v = [0:2:10 10:10:100 100:20:200];
    else
        v = [0:2:10 10:10:100 100:20:200 200:floor((max(max(F))-200)/(max(max(F))/200)):max(F) max(max(F))];
    end
    hold on
    [C,h] = contour(X1,X2,F,v,'linewidth',2);
    

    % Turns off the legend for the contour
    set( get( get( h, 'Annotation'), 'LegendInformation' ), 'IconDisplayStyle', 'off' );
    
    % Appropiate labels for the plot
    xlabel('x_1','Fontsize',14)
    ylabel('x_2','Fontsize',14)
    set(gca,'fontsize',14);
    set(gcf, 'Color', 'w');
    colorbar
    axis image
end