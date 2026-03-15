% Ejemplo de uso
% qr: tu matriz 44x2 del FEM
% qF: tu vector resultado 2x1 (ej: [1; 0.2-0.5i])
% omega: la frecuencia que obtuviste en ese punto (ej: 45.2 rad/s)

Animate_Goland_Wing(qstruct, qrlc(:, 55) + 1i*qilc(:, 55), wlc(55), 'Goland_LCO.avi');