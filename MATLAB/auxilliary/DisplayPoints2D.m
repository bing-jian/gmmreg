function [axis_limits] = DisplayPoints2D(Model, Scene, sampling, axis_limits)
%%=====================================================================
%% $RCSfile: DisplayPoints2D.m,v $
%% $Author: bing.jian $
%% $Date: 2008-11-13 16:34:29 -0500 (Thu, 13 Nov 2008) $
%% $Revision: 109 $
%%=====================================================================

set(gca,'FontSize',16,'FontName','Times','FontWeight','bold');
plot(Model(:,1),Model(:,2),'r+','MarkerSize', 8, 'LineWidth',1.5);
hold on;
plot(Scene(:,1),Scene(:,2),'bo','MarkerSize', 8, 'LineWidth',1.5);
axis equal;


if (nargin<3)
%    axis_limits = determine_border(Model, Scene);
    sampling = 0;
end

m = size(Model,1);
if (sampling>0)
    for i=1:sampling:m
        text(Model(i,1), Model(i,2), [' \leftarrow',sprintf('%d',i)]);
    end
end

m = size(Scene,1);
if (sampling>0)
    for i=1:sampling:m
        text(Scene(i,1), Scene(i,2), [' \leftarrow',sprintf('%d',i)]);
    end
end


if (nargin<4)
    axis_limits = determine_border(Model, Scene);
end
xlim(axis_limits(1,:));
ylim(axis_limits(2,:));   

pbaspect([1,1,1]);