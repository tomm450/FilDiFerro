% Load saved figures
c=hgload('03.fig');
k=hgload('04.fig');
j=hgload('05.fig');
% Prepare subplots
figure
h(1)=subplot(3,1,1);
h(2)=subplot(3,1,2);
h(3)=subplot(3,1,3);

% Paste figures on the subplots
copyobj(allchild(get(c,'CurrentAxes')),h(1));
copyobj(allchild(get(k,'CurrentAxes')),h(2));
copyobj(allchild(get(j,'CurrentAxes')),h(3));

% % Add legends
% l(1)=legend(h(1),'LegendForFirstFigure')
% l(2)=legend(h(2),'LegendForSecondFigure')

%%

close all

ax = subplot(1, 3, 1);

fig = openfig('04.fig')%, 'visible', 'off');

imh = findobj(fig, 'type', 'image');
copyobj(imh, ax);

delete(fig);

