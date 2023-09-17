figure(2);

plot([-10 L*2], [0 0], 'black', 'LineWidth',0.5);
xlim([min(x)-max(x)/10 max(x)+max(x)/10])
xlabel('Dominio')
ylim([-0.6 0.6])
ylabel('Valor de la onda')
title('Propagación de la onda');
hold on

plot([(x0-l0/2) (x0+l0/2)], [0 0], 'r', 'LineWidth',2);
hold on

plot([(x0+l0/2) (x0+l0/2)], [-1 1], 'b', 'LineWidth',1, 'LineStyle','--');
hold on

plot([(x0-l0/2) (x0-l0/2)], [-1 1], 'b', 'LineWidth',1, 'LineStyle','--');
hold on

txt = ['\uparrow x_0 = ', num2str(x0)];
text(x0,-0.05,txt)

txt = ['\leftarrow l_0 = ', num2str(l0)];
text(x0+l0/2+0.1,-0.3,txt)

txt = ['c = ', num2str(E)];
text(0,-0.2,txt)

txt = ['c = ', num2str(E)];
text(9,-0.2,txt)

txt = ['c_0 = ', num2str(E0)];
text(x0-0.4,-0.2,txt)

for n = 2 : 75 : tnts
        plot(x,U(:,n));
        hold on
end

exportgraphics(gcf,'graph2.pdf','ContentType','vector')