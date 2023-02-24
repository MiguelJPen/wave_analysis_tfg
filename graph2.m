figure(2);

plot([-10 L*2], [0 0], 'black', 'LineWidth',0.5);
xlim([min(x)-max(x)/10 max(x)+max(x)/10])
ylim([-0.6 0.6])
title('Propagación de la onda');
hold on

V = [(x0-l0) -1; (x0-l0) 1; (x0+l0) 1; (x0+l0) -1];
F = [1 2 3 4 1];
patch('Faces',F,'Vertices',V, 'FaceColor', 'r', ...
    'EdgeColor', 'none', 'FaceAlpha', 0.3);
hold on

txt = ['\leftarrow x0 = ', num2str(x0)];
text(x0,-0.8,txt)

txt = ['\leftarrow l0 = ', num2str(l0)];
text(x0+l0,-0.45,txt)

txt = ['c = ', num2str(E)];
text(0,-0.55,txt)

txt = ['c0 = ', num2str(E0)];
text(x0+0.2,-0.55,txt)

for n = 2 : 75 : tnts
        plot(x,U(:,n));
        hold on
end

exportgraphics(gcf,'graph2.pdf','ContentType','vector')