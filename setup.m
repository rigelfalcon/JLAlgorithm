% restoredefaultpath
addpath(genpath('./external/'))

addpath(genpath('./utility/'))

set(groot, 'defaultFigureCloseRequestFcn', 'delete(gcf)');
% profile -memory on;
profile('-memory','on');  % an alternative
