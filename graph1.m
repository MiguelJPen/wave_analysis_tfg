figure(1);

plot([0 L], [0 0], 'black', 'LineWidth',2);
xlim([min(x)-max(x)/10 max(x)+max(x)/10])
xlabel('Dominio')
ylim([-1 1])
ylabel('Valor de la onda')
title('Representación del medio');
hold on

plot([(x0-l0/2) (x0+l0/2)], [0 0], 'r', 'LineWidth',2);
hold on

plot([(x0+l0/2) (x0+l0/2)], [-1 1], 'b', 'LineWidth',1, 'LineStyle','--');
hold on

plot([(x0-l0/2) (x0-l0/2)], [-1 1], 'b', 'LineWidth',1, 'LineStyle','--');

txt = {'u(0,t) = f(t)', 't \in [t_0,T]'};
text(0,-0.2,txt)

txt = {'u(L,t) = g(t)', 't \in [t_0,T]'};
text(L - 1.5,-0.2,txt)

txt = {'t = t_0 + i*dt', 'i = 0, 1, ..., M', 'dt: paso de tiempo', 't_0: tiempo inicial'};
text(-0.2,0.6,txt)

txt = ['\uparrow x_0 = ', num2str(x0)];
text(x0,-0.1,txt)

txt = ['\leftarrow l_0 = ', num2str(l0)];
text(x0+l0/2+0.1,-0.5,txt)

txt = ['c = ', num2str(E)];
text(0,0.1,txt)

txt = ['c = ', num2str(E)];
text(9,0.1,txt)

txt = ['c_0 = ', num2str(E0)];
text(x0-0.4,0.1,txt)
hold off

exportgraphics(gcf,'graph1.pdf','ContentType','vector')