function createfigure1(X1, Y1, S1, C1, X2, Y2, C2, X3, Y3, C3, X4, Y4, C4, X5, Y5, C5, X6, Y6, C6, la, lb, lc, ld, le, lf)
figure('Name','Figure');
axes1 = axes('Position',[0.13 0.11 0.564273127753304 0.811854304635762]);
hold(axes1,'on');
scatter(X1,Y1,S1,C1,'Marker','.');
scatter(X2,Y2,S1,C2,'Marker','.');
scatter(X3,Y3,S1,C3,'Marker','.');
scatter(X4,Y4,S1,C4,'Marker','.');
scatter(X5,Y5,S1,C5,'Marker','.');
scatter(X6,Y6,S1,C6,'Marker','.');
xlabel('Bioactivity pORACL A','FontName','Arial');
ylabel('Bioactivity pORACL B','FontName','Arial');
xlim(axes1,[1e-16 1]);
ylim(axes1,[1e-16 1]);
legend([la lb lc ld le lf],'Location','southoutside');
box(axes1,'on');
% Set the remaining axes properties
set(axes1,'PlotBoxAspectRatio',[1 1 1],'Box','on','FontName','Arial','FontSmoothing','off','XGrid','on',...
    'XMinorTick','on','XScale','log','XTick',[1e-16 1e-06 0.001 1],'YGrid','on',...
    'YMinorTick','on','YScale','log','YTick',[1e-16 1e-06 0.001 1],'GridColor','k','GridAlpha',0.5);
