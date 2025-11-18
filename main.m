%% ========== CÓDIGO PRINCIPAL ==========

clear; close all; clc;

% =========== Carga de Parámetros y Perturbaciones ==========
p = parametros();
pert = Entradas();

tspan = 0:69:24*3600;
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
test_temp_amb = @(t) 30;
test_irr_sol = @(t) 800;
test_vel_vient = @(t) 5.0;

% --- Fase 1: Estabilizar en 45°C ---
p_test.ref = 45 + 273.15;
test_y0_1 = [30 + 273.15, 30 + 273.15];
tspan_1 = 0:60:1*3600;

odefun_test_1 = @(t, y) modelo(t, y, p_test, test_temp_amb, ...
                                 test_irr_sol, test_vel_vient, ...
                                 1, @(t) 0, @(t) 0, test_Kp);
[t_test_1, y_test_1] = ode45(odefun_test_1, tspan_1, test_y0_1);
T_p_test_1 = y_test_1(:, 1) - 273.15;

% --- CÁLCULO DE VOLTAJE (FASE 1) ---
% (Debe hacerse ANTES de cambiar p_test.ref)
salida_sin_limite_1 = test_Kp * (p_test.ref - y_test_1(:,1)) + p_test.offset;
salida_vent_1 = min(p_test.V_vent_MAX, max(salida_sin_limite_1, p_test.V_MIN));
salida_bomb_1 = min(p_test.V_bomb_suficiente, max(salida_sin_limite_1, p_test.V_MIN));
 
% --- Fase 2: Aplicar Escalón a 55°C ---
p_test.ref = 55 + 273.15;
tspan_2 = 1*3600:60:2*3600;
test_y0_2 = y_test_1(end, :);

odefun_test_2 = @(t, y) modelo(t, y, p_test, test_temp_amb, ...
                                 test_irr_sol, test_vel_vient, ...
                                 1, @(t) 0, @(t) 0, test_Kp);
[t_test_2, y_test_2] = ode45(odefun_test_2, tspan_2, test_y0_2);
T_p_test_2 = y_test_2(:, 1) - 273.15;

% --- CÁLCULO DE VOLTAJE (FASE 2) ---
salida_sin_limite_2 = test_Kp * (p_test.ref - y_test_2(:,1)) + p_test.offset;
salida_vent_2 = min(p_test.V_vent_MAX, max(salida_sin_limite_2, p_test.V_MIN));
salida_bomb_2 = min(p_test.V_bomb_suficiente, max(salida_sin_limite_2, p_test.V_MIN));

% --- CÁLCULO DE MÉTRICAS (OBJETIVO 5) ---
fprintf('--- Métricas del Controlador (Step Test, Kp=%.2f) ---\n', test_Kp);

% Definir las referencias del escalón
Ref_C_inicial = 40; % El valor de p_test.ref en Fase 1
Ref_C_final = 55;   % El valor de p_test.ref en Fase 2
T_p_fase2_C = y_test_2(:,1) - 273.15; % Vector de T de la Fase 2

% 1. Error de Estado Estacionario (ESS)
%    (El hallazgo más importante)
T_estable_real = T_p_fase2_C(end); % Temperatura final real (tus ~48°C)
ESS_C = Ref_C_final - T_estable_real;
fprintf('Error de Estado Estacionario (ESS): %.2f °C\n', ESS_C);

% 2. Sobrepaso (Overshoot)
T_max = max(T_p_fase2_C);
% El sobrepaso es relativo al valor ESTABLE, no a la referencia
Sobrepaso_C = T_max - T_estable_real; 
fprintf('Sobrepaso (relativo a %.1f°C): %.2f °C\n', T_estable_real, Sobrepaso_C);

% 3. Tiempo de Estabilización (Settling Time)
%    (Calculado alrededor del NUEVO punto estable: 48°C)

% Definimos la banda de tolerancia (ej. 2%) con respecto al valor estable
banda_tolerancia = 0.02 * T_estable_real; % 2% de 48°C ≈ 0.96°C
banda_superior = T_estable_real + banda_tolerancia;
banda_inferior = T_estable_real - banda_tolerancia;

% Encontrar el *último* índice (fila) donde la Temp ESTÁ FUERA de la banda
indice_fuera = find(T_p_fase2_C > banda_superior | T_p_fase2_C < banda_inferior, 1, 'last');
if isempty(indice_fuera) % Si nunca salió de la banda
    indice_fuera = 1; 
end
T_estabilizacion_seg = t_test_2(indice_fuera) - tspan_1(end); 
fprintf('Tiempo de Estabilización (banda 2%% alrededor de %.1f°C): %.1f segundos\n', T_estable_real, T_estabilizacion_seg);

% 4. Tiempo de Subida (Rise Time)
%    (Esto sigue igual, 10% a 90% del viaje de 40°C a 48°C)
Amplitud_Viaje_Real = T_estable_real - Ref_C_inicial;
T_10_pct = Ref_C_inicial + (0.1 * Amplitud_Viaje_Real);
T_90_pct = Ref_C_inicial + (0.9 * Amplitud_Viaje_Real);
indice_t10 = find(T_p_fase2_C >= T_10_pct, 1, 'first');
indice_t90 = find(T_p_fase2_C >= T_90_pct, 1, 'first');
T_subida_seg = t_test_2(indice_t90) - t_test_2(indice_t10);
fprintf('Tiempo de Subida (10-90%% del viaje real): %.1f segundos\n', T_subida_seg);

% --- Gráfico del Step Test (Actualizado con Métricas) ---
figure("Name", "Step Test (Objetivo 5)", "NumberTitle", "off");

% Gráfico de Temperatura
subplot(1, 2 , 1);
plot(t_test_1/3600, T_p_test_1, "LineWidth", 2, "Color", [0.5 0.5 0.5]); % Fase 1 en gris
hold on; grid on;
plot(t_test_2/3600, T_p_test_2, "LineWidth", 2, "Color", "b"); % Fase 2 en azul
% --- Líneas de Métricas ---
yline(Ref_C_final, 'r--', 'Referencia Final');
yline(T_max, 'k:', sprintf('Sobrepaso (%.2f°C)', T_max));
xline(t_test_2(indice_fuera)/3600, 'g:', sprintf('T. Estab. (%.1fs)', T_estabilizacion_seg));
% ---
title("Respuesta al Escalón (40°C a 55°C)");
legend('Pre-estabilización', 'Respuesta al Escalón', 'Location', 'southeast');
xlim([0, 2]); xticks(0:0.5:2);
ylim([35, 60]); % Ajustar Ylim
xlabel("Tiempo (h)"); ylabel("Temperatura (°C)");

% Gráfico de Voltaje
subplot(1, 2, 2);
plot(t_test_1/3600, salida_vent_1, "LineStyle", "--", "LineWidth", 2, "Color", [0.5 0.5 1]); % Azul claro
hold on; grid on;
plot(t_test_2/3600, salida_vent_2, "LineStyle", "--", "LineWidth", 2, "Color", "b"); % Azul fuerte
plot(t_test_1/3600, salida_bomb_1, "LineStyle", "-.", "LineWidth", 2, "Color", [1 0.5 0.5]); % Rojo claro
plot(t_test_2/3600, salida_bomb_2, "LineStyle", "-.", "LineWidth", 2, "Color", "r"); % Rojo fuerte
title("Voltaje del Controlador");
xlim([0,2]); xticks(0:0.5:2);
xline(1, 'g--', 'Cambio de Ref.');
ylim([-1, 13]);
xlabel("Tiempo (h)"); ylabel("Voltaje (V)");

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