function graficos(nombre, t, T_p, temp_amb, irr_sol, sal_vent, ...
                    sal_bomb)

    % ---------- Temperatura Panel Solar ----------
    figure("Name", nombre, "NumberTitle", "off");
    subplot(2, 2, 1);
    plot(t/3600, T_p);
    title("Temperatura del Panel Solar");
    xlim([0,24]); xticks(0:2:24);
    ylim([0, 60]);
    xlabel("Tiempo (h)"); ylabel("Temperatura (°C)");
    grid on;
    
    % ---------- Voltajes Ventilador/Bomba ----------
    subplot(2, 2, 2);
    plot(t/3600, sal_vent, "LineStyle", "--", "LineWidth", 2, "Color", "b");
    hold on; grid on;
    plot(t/3600, sal_bomb, "LineStyle", "-.", "LineWidth", 2, "Color", "r");
    legend("V. Ventilador", "V. Bomba", "Location", "northwest");
    title("Voltaje aplicado a la bomba y ventilador");
    xlim([0,24]); xticks(0:2:24);
    ylim([-1, 13]);
    xlabel("Tiempo (h)"); ylabel("Voltaje (V)");

    % ---------- Temperatura Ambiente ----------

    subplot(2, 2, 3);
    plot(t/3600, temp_amb, "LineWidth", 1, "Color", "g");
    title("Temperatura Ambiente");
    xlim([0,24]); xticks(0:2:24);
    ylim([0, 40]);
    xlabel("Tiempo (h)"); ylabel("Temperatura (°C)");
    grid on;

    % ---------- Irradiancia Solar ----------
    subplot(2, 2, 4);
    plot(t/3600, irr_sol, "LineWidth", 1, "Color", "y");
    title("Irradiancia Solar");
    xlim([0,24]); xticks(0:2:24);
    ylim([-100, 1100]);
    xlabel("Tiempo (h)"); ylabel("Irradiancia (W/m^2)");
    grid on;

end