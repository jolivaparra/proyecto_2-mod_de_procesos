%% ========== GRÁFICAS DE LOS PERFILES ==========
function graficos_perfil(temp_amb, V_bomb, V_vent, t_ctrl, T_p_ctrl, eta_ctrl, t, T_p, T_w_in, eta, i)

    figure('Name', ['Salidas - Perfil ',num2str(i)], 'NumberTitle','off');
    % GRAFICA: Panel vs Ambiente vs Panel s/n refrigeración
    subplot(2,2,1);
    plot(t_ctrl/3600, T_p_ctrl -273.15, '-.m', 'LineWidth', 1);  hold on;
    plot(t/3600, temp_amb, "--", "LineWidth", 1, 'Color', [0 0.5 0]);
    plot(t/3600, T_p - 273.15, 'r', "LineWidth", 2);

    xlabel("Tiempo (h)"); ylabel("Temperatura (°C)");
    title("(A) Temperatura del Panel vs Ambiente");
    legend("s/n Refr.", "Amb", "c/n Refr.", "Location", "northwest");
    xticks(0:2:24);
    xlim([0,24]);
    grid on;

    % GRAFICA: Cambio en eficiencia
    subplot(2, 2, 2);
    plot(t_ctrl/3600, eta_ctrl*100, "--", "LineWidth", 1, 'Color', [0 0.5 0]); hold on;
    plot(t/3600, eta*100, "-b", "LineWidth", 2);
    xlabel("Tiempo (h)"); ylabel("Eficiencia (%)");
    title("(B) Eficiencia del Panel Solar");
    legend("s/n Refr. Activa", "c/n Refr. Activa", 'Location', 'southeast')
    xticks(0:2:24);
    xlim([0,24]);
    grid on;

    % GRAFICA: Agua
    subplot(2,2,3);
    plot(t/3600, T_w_in - 273.15, "-b", "LineWidth", 2);
    xlabel("Tiempo (h)"); ylabel("Temperatura (°C)");
    title("(C) Temperatura del Agua entrante");
    xticks(0:2:24);
    xlim([0,24]);
    grid on;

    % GRAFICA: Voltajes
    subplot(2,2,4);
    plot(t/3600, V_bomb, "-.b", "LineWidth", 2); hold on;
    plot(t/3600, V_vent, "--r", "LineWidth", 2);
    xlabel("Tiempo (h)"); ylabel("Voltaje (V)");
    title("(D) Voltajes de Entrada");
    legend("Bomba", "Ventilador", "Location", "northwest");
    xticks(0:2:24); yticks(0:6:25); 
    xlim([0,24]); ylim([-1 26]);
    grid on;

end