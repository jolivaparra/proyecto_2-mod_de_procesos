%% ========== CÓDIGO PRINCIPAL ==========

clear; close all; clc;

% =========== Carga de Parámetros ==========
p = parametros();

%% ========== VARIABLES DE ENTRADA ==========

% Perfil de Entrada: Escalón

voltaje_bomba_escalon = @(t) 24 * ((9 <= t/3600) & (t/3600 <= 19));
voltaje_ventilador_escalon = @(t) 12 * ((9 <= t/3600) & (t/3600 <= 19));

% Perfil de entrada: Rampa

voltaje_bomba_rampa = @(t) 1 * (t/3600);
voltaje_ventilador_rampa = @(t) (1/2) * (t/3600);

% Perfil de entrada: Sinusoidal

voltaje_bomba_sinusoidal = @(t) 12 + 12 * cos(pi/2 * (t/3600)-pi);
voltaje_ventilador_sinusoidal = @(t) 6 + 6 * cos(pi/2 * (t/3600)-pi);

%% ========== VARIABLES EXTERNAS ==========

% Temperatura ambiente
function T_amb = temperatura_ambiente(t)
   
    t_h = t/3600;

    t_min = 6.5;
    t_sunset = 20.5;

    Y = t_sunset-t_min;
    Z = 24-Y;

    T_x = 35.2;
    T_N = 10;
    T_s = 25.1018;

    a = 1.8;
    b = 2.5;

    if t_h>=t_min && t_h<=t_sunset
        
        T_amb = (T_x-T_N)*sin(pi*(t_h-t_min)/(Y+2*a))+T_N;
   
    else

        % Ajuste de continuidad y periodicidad 
        % para comenzar simulaciones en tiempo 00:00
        if t_h<t_min
            T_amb = T_N+(T_s-T_N)*exp(-(b*((t_h+24)-t_sunset))/(Z));
        
        else
            T_amb = T_N+(T_s-T_N)*exp(-(b*(t_h-t_sunset))/(Z));
        end
    end
end

% Irradiancia solar
function G = irradiancia_solar(t_segundos)

    G_max = 1000;      % Irradiancia pico (W/m^2)
    t_amanecer = 6.5;    % Hora del amanecer (horas)
    t_atardecer = 20.5;  % Hora del ocaso (horas)
    t_horas = t_segundos / 3600;

    arg = pi * (t_horas-t_amanecer)/(t_atardecer-t_amanecer);

    G = (t_amanecer<=t_horas & t_horas<=t_atardecer).*G_max.*sin(arg);

end

% Velocidad del viento
function v_vient = velocidad_viento(t)
    v_vient = 5+3*sin(2*pi*t/3600*24);
end

%% ========== MODELO ==========
function dydt = modelo(t, y, p, V_bomb_fun, V_vent_fun)

    % Variables de estado en [K]
    T_p = y(1);
    T_w_in = y(2);

    % Variables de Entrada
    V_bomb = V_bomb_fun(t);
    V_vent = V_vent_fun(t);

    % Variables Externas
    v_vient = velocidad_viento(t);
    T_amb = temperatura_ambiente(t) + 273.15;
    G = irradiancia_solar(t);
    
    % Temperatura del cielo
    T_cielo = T_amb - 8;
    
    % Flujo másico
    m_dot_w = p.rho_w * p.k_Qf * V_bomb;

    % Convección panel solar
    if v_vient<5
        h_p = 5.7+3.8*v_vient;
    else
        h_p = 6.47+v_vient^(0.78);
    end
    

    if V_bomb > 0
        % ---------- BOMBA ENCENDIDA ----------

        T_w_out = (2*m_dot_w * p.C_p_w * T_w_in + p.h_conv_pw * ...
                  p.A_c * (2*T_p - T_w_in)) / (2*m_dot_w * p.C_p_w ...
                  + p.h_conv_pw * p.A_c);
        
        % DeltaT_ml sin indeterminaciones
        % Para evitar posibles problemas numéricos
        if abs(T_w_out - T_w_in) > 0 && (T_w_out - T_amb) > 0 && ...
                (T_w_in - T_amb) > 0

            arg_log = (T_w_out - T_amb) / (T_w_in - T_amb);

            if arg_log > 0 && abs(arg_log - 1) > 0
                Delta_T_ml = (T_w_out - T_w_in) / log(arg_log);

            else
                Delta_T_ml = T_w_in - T_amb;
            end
        else
            Delta_T_ml = T_w_in - T_amb;
        end
        
    else
        % ---------- BOMBA APAGADA ----------
        T_w_out = T_w_in;  % Sin flujo, no hay gradiente
        Delta_T_ml = T_w_in - T_amb;  % Para agua estancada, diferencia simple
    end
    

    % Flujos de calor
    Q_sol = p.A_p * G * (p.alpha_p - p.eta_ref * ...
             (1 + (p.mu_Voc/p.V_mp)*(T_p - p.T_ref)));

    Q_conv_p = h_p * p.A_p * (T_p - T_amb);

    Q_rad_p = p.epsilon_p * p.sigma * p.A_p * (T_p^4 - T_cielo^4);

    Q_cool = m_dot_w * p.C_p_w * (T_w_out - T_w_in);

    Q_int = (p.U_base + p.k_vent * V_vent) * p.A_int * Delta_T_ml;

    % Ecuaciones diferenciales
    dT_p_dt = (Q_sol - Q_conv_p - Q_rad_p - Q_cool) / (p.m_p * p.C_p_p);
    dT_w_in_dt = (Q_cool - Q_int) / (p.m_w * p.C_p_w);

    % Sistema
    dydt = [dT_p_dt;
            dT_w_in_dt];
end

%% ========== SIMULACIÓN Y GRÁFICOS ==========

% ---- Condiciones iniciales -----
T_amb_0 = temperatura_ambiente(0) + 273.15; % [K]
y0 = [T_amb_0, T_amb_0];  % Panel ligeramente más caliente
tspan = 0:24*3600;

% Temperatura del panel sin refrigeración para comparar
tspan_ctrl = 0:450:24*3600;
[t_ctrl, y_ctrl] = ode45(@(t, y) modelo(t, y, p, @(t) 0, @(t) 0), tspan_ctrl, y0);
T_p_ctrl = y_ctrl(:, 1);
T_w_in_ctrl = y_ctrl(:, 2);
eta_ctrl = arrayfun(@(T_p) p.eta_ref * (1 + (p.mu_Voc/p.V_mp)*(T_p - p.T_ref)), T_p_ctrl); 

% TABLA: |Temp. Máx|Hora Temp. Max (Hr.Mn)|Temp. Media durante día|
%                               |Reducción de Temp|Aumento de eff. panel|

% Temperatura y momento cuando la Temp está en su momento máx
[T_max_ctrl, idx_T_max_ctrl] = max(T_p_ctrl);

% Se obtiene la hora exacta en formato Hr:Mn
hora_T_max_ctrl = t_ctrl(idx_T_max_ctrl)/3600;

h_T_max_ctrl = fix(hora_T_max_ctrl); %Horas
m_T_max_ctrl = mod(hora_T_max_ctrl, 1)*(60/100);%Minutos
hora_tabla_T_max_ctrl = h_T_max_ctrl + m_T_max_ctrl; % Unión

% Temperaturas medias de temperatura y eficiencia durante horas cálidas
media_T_ctrl_dia = mean(T_p_ctrl(10.5*3600<=t_ctrl & t_ctrl<=20.5*3600)) - 273.15;
media_eta_ctrl_dia = mean(eta_ctrl(10.5*3600<=t_ctrl & t_ctrl<=20.5*3600));

% Tabl con | Temp. Máx | Hora Temp. Max (Hr.Mn) | Temp. Media durante día |
datos_perfiles = table(T_max_ctrl-273.15, hora_T_max_ctrl, media_T_ctrl_dia, ...
                        0, 0, ...
                       'VariableNames', {'T_max', 'hora_T_max', ...
                       'media_T_dia', 'reduccion_T', 'reduccion_eta'});


% Perfiles de entrada al sistema
perfiles_entrada = {
    voltaje_bomba_escalon,      voltaje_ventilador_escalon; 
    voltaje_bomba_rampa,        voltaje_ventilador_rampa;
    voltaje_bomba_sinusoidal,   voltaje_ventilador_sinusoidal
};

% Generación de gráficos por cada perfil de entrada
for i = 1:size(perfiles_entrada,1)
    
    % Se determina cada variable de entrada
    V_bomb_fun = perfiles_entrada{i, 1};
    V_vent_fun = perfiles_entrada{i, 2};

    % Resolución del sistema
    [t, y] = ode45(@(t, y) modelo(t, y, p, V_bomb_fun, V_vent_fun), tspan, y0);
    T_p = y(:, 1);
    T_w_in = y(:, 2);

    % Para graficar variables de entrada y externas
    temp_amb    = arrayfun(@(t) temperatura_ambiente(t), t);
    V_bomb      = arrayfun(@(t) V_bomb_fun(t), t);
    V_vent      = arrayfun(@(t) V_vent_fun(t), t);
    eta         = arrayfun(@(T_p) p.eta_ref * (1 + (p.mu_Voc/p.V_mp)*(T_p - p.T_ref)), T_p);

    % Gráficos
    graficos_perfil(temp_amb, V_bomb, V_vent, t_ctrl, T_p_ctrl, eta_ctrl, t, T_p, T_w_in, eta, i);


    % ----------- TABLA -----------
    % Añadir datos relevantes a la tabla
    [T_max, idx_T_max] = max(T_p); % Temp. y momento en T. Máx
    t_h_T_max = t(idx_T_max)/3600; % Tiempo a horas

    horas_T_max = fix(t_h_T_max); % Hora
    min_T_max = mod(t_h_T_max, 1)*(60/100); % Minutos
    hora_T_max = horas_T_max + min_T_max; % Unión
    
    % Medias de Temp. y Eficiencia durante horas de mayor calor
    media_T_dia = mean(T_p(10.5*3600<=t & t<=20.5*3600)) - 273.15;
    media_eta_dia = mean(eta(10.5*3600<=t & t<=20.5*3600));
    
    % Cálculo de % de reducción y aumento
    reduccion_T = (media_T_ctrl_dia - media_T_dia)/media_T_ctrl_dia * 100;
    aumento_eta = (media_eta_dia - media_eta_ctrl_dia)/media_eta_ctrl_dia * 100;
    
    datos_perfiles(end+1, :) = {T_max-273.15, hora_T_max, media_T_dia, ...
                                reduccion_T, aumento_eta};
    % ------------------------------
end

%% GRAFICOS DE ENTRADAS EXTERNAS

tspan = 0:24*3600;

% Para graficar
T_amb_valores = arrayfun(@temperatura_ambiente, tspan);
G_valores = arrayfun(@irradiancia_solar, tspan);

% Graficos de entradas externas
figure('Name', 'Temperatura Ambiente', 'NumberTitle','off');
subplot(1, 2, 1);
plot(tspan / 3600, T_amb_valores, 'LineWidth', 2);
xlabel('Tiempo (h)');
ylabel('Temperatura (°C)');
title('Temperatura ambiente durante el día');
xticks(0:2:24); yticks(0:5:40);
xlim([0,24]); ylim([7.5, 40]);
grid on;

figure('Name', 'Irradiancia Solar', 'NumberTitle','off');
subplot(2, 1, 1);
plot(tspan / 3600, G_valores, "LineWidth",2);
xlabel('Tiempo (h)');
ylabel('Irradiancia solar (W/m^2)');
title('Irradiancia solar durante un día despejado');
xticks(0:2:24);
xlim([0,24]); ylim([-100, 1050]);
grid on
