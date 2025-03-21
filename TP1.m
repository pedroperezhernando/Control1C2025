optionss=bodeoptions;
optionss.MagVisible='off';
optionss.PhaseMatching='on';
optionss.PhaseMatchingValue=-180;
optionss.PhaseMatchingFreq=1;
optionss.Grid='on';
optionss.MagVisible='on';
%%
%Puntos de modelado
p_trab = [0.10,0.20,0.30,0.40,0.50,0.60,0.70,0.80]

% Inicializar un array de celdas para almacenar las funciones de transferencia
res = cell(size(p_trab));

% Loop sobre cada punto de trabajo
for i = 1:length(p_trab)
    res{i} = planta(p_trab(i)); % Llamar a la función correctamente
end
   
figure; % Crear una nueva figura
hold on; % Mantener el mismo gráfico para múltiples bode

for i = 1:length(res)
    bode(res{i}); % Graficar cada función de transferencia
end

hold off; % Liberar la retención del gráfico
grid on; % Activar la cuadrícula
legend(arrayfun(@(h) sprintf('h_0 = %.2f', h), p_trab, 'UniformOutput', false)); % Agregar leyenda
title('Comparación de Diagramas de Bode para Distintos h_0');
%%

s=tf('s')

P_dis = planta(0.45)

C_PI = -(s+0.1)/(s)

L = P_dis * C_PI

bode(L)
k = db2mag(30.5)

C_PI_Prima = C_PI * k

L_prima = P_dis * C_PI_Prima
figure
bode(L_prima)
T = feedback(L_prima,0.45)
Ts = 0.014*(5)/(2*pi)
C_dig = c2d(C_PI_Prima,Ts,'tustin')

t = 0:0.01:500; % Tiempo de simulación de 0 a 10 segundos con paso de 0.01s
[y, t] = step(0.6 * T, t); % Respuesta del sistema a un escalón de 0.6

figure;
plot(t, y, 'r', 'LineWidth', 2);
xlabel('Tiempo (s)');
ylabel('Salida');
title('Respuesta al Escalón de Amplitud 0.6');
grid on;


function [P1] = planta(h_0)
    x1_e = h_0
    q_i = 8/60 % Litros/seg caudal de entrada cte
    d_2 = 10.65/1000 %metros diametro canio de salida
    l_1 = 10/100%metros lado base inferior
    l_2 = 40/100%metros lado base superior
    l = 0.9%metros de la base del tanque al tope.
    g = 9.81%metros/seg**2, gravedad Asumo que x1_punto eq = 0 (no hay variacion de altura)
    u1_e = (q_i)/(pi*(d_2^2/4)*sqrt(2*g*x1_e))
    u1_e = u1_e/1000 % Conversion de dm^3/m^3
    y1_e = x1_e;
    s = tf('s');
    syms x1 u y;

    x = x1;
    f = (q_i/1000 - (u * pi * (d_2^2 / 4) * sqrt(2 * g * x1))) / (l_1 + (l_2 - l_1) * (x1 / l) )^2;
    y = x1 %Paso a decimetro ^3 de m^3 %Litros 

    A = jacobian(f,x);
    B = jacobian(f,u);
    C = jacobian(y,x);
    D = jacobian(y,u);

    A = subs(A, str2sym({'x1','u','y'}), {x1_e,u1_e,y1_e});
    B = subs(B, str2sym({'x1','u','y'}), {x1_e,u1_e,y1_e});
    C = subs(C, str2sym({'x1','u','y'}), {x1_e,u1_e,y1_e});
    D = subs(D, str2sym({'x1','u','y'}), {x1_e,u1_e,y1_e});

    Ass = double(A);
    Bss = double(B);
    Css = double(C);
    Dss = double(D);

    % Transferencia
    P1 = zpk(ss(Ass, Bss, Css, Dss));
end