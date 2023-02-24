figure(1);

plot([0 L], [0 0], 'black', 'LineWidth',2);
xlim([min(x)-max(x)/10 max(x)+max(x)/10])
ylim([-1 1])
hold on

plot([x0 x0], [-1 1], 'r', 'LineWidth',1);
hold on

Vr = [(x0-l0) -1; (x0-l0) 1; (x0+l0) 1; (x0+l0) -1];
F = [1 2 3 4 1];
patch('Faces',F,'Vertices',Vr, 'FaceColor', 'r', ...
    'EdgeColor', 'none', 'FaceAlpha', 0.3);

txt = {'u(0,t) = f(t)', 't \in [t_0,T]'};
text(0,-0.2,txt)

txt = {'u(L,t) = g(t)', 't \in [t_0,T]'};
text(L - 1.5,-0.2,txt)

txt = {'t = t_0 + i*dt', 'i = 0, 1, ..., M', 'dt: paso de tiempo', 't_0: tiempo inicial'};
text(0,0.6,txt)

txt = ['\leftarrow x0 = ', num2str(x0)];
text(x0,-0.8,txt)

txt = ['\leftarrow l0 = ', num2str(l0)];
text(x0+l0,-0.65,txt)

txt = ['c = ', num2str(E)];
text(0,-0.9,txt)

txt = ['c0 = ', num2str(E0)];
text(x0+0.2,-0.9,txt)
hold off

exportgraphics(gcf,'graph1.pdf','ContentType','vector')