figure(3);

plot((0.1 : 0.1 : 10), result(1, :), 'b');
xlim([0 10.5])
xlabel('Tiempo')
ylim([min(result)-max(result)/10 max(result)+max(result)/10])
ylabel('Valor de la onda')
title('u(L, t) = g(t)');
hold on

exportgraphics(gcf,'graph3.pdf','ContentType','vector')