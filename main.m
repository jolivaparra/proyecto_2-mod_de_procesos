%% ========== CÓDIGO PRINCIPAL ==========

clear; close all; clc;

% =========== Carga de Parámetros y Perturbaciones ==========
p = parametros();
pert = Entradas();

tspan = 0:24*3600;
T0 = pert.temperatura_ambiente_despejado(0) + 273.15;
y0 = [T0, T0];

%% ========== SIN REFRIGERACION ==========

odefun_snrefr = @(t, y) modelo(t, y, p, ...
                                @(t) pert.temperatura_ambiente_despejado(t), ...
                                @(t) pert.irradiancia_solar_despejado(t), ...
                                @(t) pert.velocidad_viento(t), ...
                                0, ...
                                @(t) 0, ...
                                @(t) 0, ...
                                0);

[t_sr, y_sr] = ode45(odefun_snrefr, tspan, y0);
T_p_sr = y_sr(:, 1) - 273.15;

figure("Name", "Sin Refrigeración", "NumberTitle", "off");
plot(t_sr/3600, T_p_sr, "LineWidth", 2, "Color", "b");
title("Temperatura del Panel, sin uso de refrigeración");
xlim([0, 24]); xticks(0:2:24);
grid on;

%% ========== LAZO ABIERTO (ESCALÓN)  ==========

voltaje_bomba_escalon = @(t) 24 * ((9 <= t/3600) & (t/3600 <= 19));
voltaje_ventilador_escalon = @(t) 12 * ((9 <= t/3600) & (t/3600 <= 19));

odefun_la = @(t, y) modelo(t, y, p, ...
                                @(t) pert.temperatura_ambiente_despejado(t), ...
                                @(t) pert.irradiancia_solar_despejado(t), ...
                                @(t) pert.velocidad_viento(t), ...
                                0, ...
                                voltaje_bomba_escalon, ...
                                voltaje_ventilador_escalon, ...
                                0);

[t_la, y_la] = ode45(odefun_la, tspan, y0);
T_p_la = y_la(:, 1) - 273.15;

figure("Name", "Lazo Abierto", "NumberTitle", "off");
plot(t_la/3600, T_p_la, "LineWidth", 2, "Color", "b");
title("Temperatura del Panel: Control en Lazo Abierto (Escalón)");
xlim([0, 24]); xticks(0:2:24);
ylim([0, 60]);
grid on;

%% ========== TEST EN ENTORNO CONTROLADO (PARA ANÁLISIS DE ESTABILIDAD) ==========
p_test = p;
test_Kp = p_test.K_p(1);
test_temp_amb = @(t) 25;
test_irr_sol = @(t) 800;
test_vel_vient = @(t) 5.0;
p_test.ref = 40 + 273.15;

test_y0_1 = [30 + 273.15, 30 + 273.15];
tspan_1 = 0:1*3600;

odefun_test_1 = @(t, y) modelo(t, y, p_test, test_temp_amb, ...
                                test_irr_sol, ...
                                test_vel_vient, ...
                                1, ...
                                @(t) 0, ...
                                @(t) 0, ...
                                test_Kp);

[t_test_1, y_test_1h] = ode45(odefun_test_1, tspan_1, test_y0_1);
T_p_test_1 = y_test_1h(:, 1) - 273.15;
                            
p_test.ref = 55 + 273.15;
tspan_2 = 1*3600:2*3600;
test_y0_2 = y_test_1h(end, :);
p.ref = 45 + 273.15;

odefun_test_2 = @(t, y) modelo(t, y, p_test, test_temp_amb, ...
                                test_irr_sol, ...
                                test_vel_vient, ...
                                1, ...
                                @(t) 0, ...
                                @(t) 0, ...
                                test_Kp);

[t_test_2, y_test_2] = ode45(odefun_test_2, tspan_2, test_y0_2);
T_p_test_2h = y_test_2(:, 1) - 273.15;

figure("Name", "Step Test", "NumberTitle", "off");
plot(t_test_1/3600, T_p_test_1, "LineWidth", 2, "Color", "b");
hold on; grid on;
plot(t_test_2/3600, T_p_test_2h);
title("Temperatura del Panel en entorno controlado.");
ylim([0, 60]);
grid on;

%% ========== LAZO CERRADO ==========
nombre_perfiles = ["Despejado", "Nublado", "Parcialmente Nublado"];

perfil_pert = {
    @(t) pert.temperatura_ambiente_despejado(t), @(t) pert.irradiancia_solar_despejado(t);
    @(t) pert.temperatura_ambiente_nublado(t), @(t) pert.irradiancia_solar_nublado(t);
    @(t) pert.temperatura_ambiente_intermitente(t), @(t) pert.irradiancia_solar_intermitente(t)
};

for i=1:size(perfil_pert, 1)
    
    temperatura_ambiente = perfil_pert{i, 1};
    irradiancia_solar = perfil_pert{i, 2};
    velocidad_viento = @(t) pert.velocidad_viento(t);

    if i==1
        for j=1:size(p.K_p , 2)

            nombre_graf = "Día Despejado: Kp "+ num2str(j);
            
            graficos(nombre_graf, p, temperatura_ambiente, irradiancia_solar, ...
                    velocidad_viento, 1, @(t) 0, @(t) 0, p.K_p(j));
        
        end

    else
        
        graficos(nombre_perfiles(i), p, temperatura_ambiente, irradiancia_solar, ...
                    velocidad_viento, 1, @(t) 0, @(t) 0, p.K_p(1));

    end
    
end