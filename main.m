%% ========== CÓDIGO PRINCIPAL ==========

clear; close all; clc;

% =========== Carga de Parámetros y Perturbaciones ==========
p = parametros();
pert = Entradas();

% ========== ENTRADA EN LAZO ABIERTO (ESCALÓN) ==========

voltaje_bomba_escalon = @(t) 24 * ((9 <= t/3600) & (t/3600 <= 19));
voltaje_ventilador_escalon = @(t) 12 * ((9 <= t/3600) & (t/3600 <= 19));

%% ========== SIMULACIONES ==========

tspan = 0:24*3600;
T0 = pert.temperatura_ambiente_despejado(0) + 273.15;
y0 = [T0, T0];

% ========== SIN REFRIGERACION ==========

odefun_sin_refrigeracion = @(t, y) modelo(t, y, p, ...
                                @(t) pert.temperatura_ambiente_despejado(t), ...
                                @(t) pert.irradiancia_solar_despejado(t), ...
                                @(t) pert.velocidad_viento(t), ...
                                0, ...
                                @(t) 0, ...
                                @(t) 0, ...
                                0);

[t_sin_refrigeracion, y_sin_refrigeracion] = ode45(odefun_sin_refrigeracion, tspan, y0);
T_p_sin_refrigeracion = y_sin_refrigeracion(:, 1) - 273.15;

figure;
plot(t_sin_refrigeracion/3600, T_p_sin_refrigeracion);
grid on;

% ========== LAZO ABIERTO ==========

odefun_lazo_abierto = @(t, y) modelo(t, y, p, ...
                                @(t) pert.temperatura_ambiente_despejado(t), ...
                                @(t) pert.irradiancia_solar_despejado(t), ...
                                @(t) pert.velocidad_viento(t), ...
                                0, ...
                                voltaje_bomba_escalon, ...
                                voltaje_ventilador_escalon, ...
                                0);

[t_lazo_abierto, y_lazo_abierto] = ode45(odefun_lazo_abierto, tspan, y0);
T_p_lazo_abierto = y_lazo_abierto(:, 1) - 273.15;

figure;
plot(t_lazo_abierto/3600, T_p_lazo_abierto);
xlim([0,24]); xticks(0:2:24);
xlabel("Tiempo (h)"); ylabel("Temperatura (°C)");
grid on;


% ========== LAZO CERRADO ==========
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