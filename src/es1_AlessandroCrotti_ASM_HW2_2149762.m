% Definizione delle coordinate nodali
nodes = [
    0, 0, 0;
    3, 0, 0;
    4, 5, 0;
    1, 5, 0;
    0, 0, 4;
    3, 0, 4;
    4, 5, 4;
    1, 5, 4
];

% Proprietà del materiale (esempio: acciaio)
E = 210e3;  % Modulo di Young in MPa
nu = 0.3;   % Poisson's ratio

% Matrice costitutiva D per materiale isotropo lineare elastico
D = (E / ((1 + nu) * (1 - 2 * nu))) * [
    1 - nu,     nu,     nu,       0,          0,          0;
        nu, 1 - nu,     nu,       0,          0,          0;
        nu,     nu, 1 - nu,       0,          0,          0;
         0,      0,      0, (1 - 2*nu)/2,     0,          0;
         0,      0,      0,       0, (1 - 2*nu)/2,      0;
         0,      0,      0,       0,          0, (1 - 2*nu)/2
];

% Punti di Gauss (2x2x2 quadrature)
gauss_pts = [
    -1/sqrt(3), -1/sqrt(3), -1/sqrt(3);
     1/sqrt(3), -1/sqrt(3), -1/sqrt(3);
     1/sqrt(3),  1/sqrt(3), -1/sqrt(3);
    -1/sqrt(3),  1/sqrt(3), -1/sqrt(3);
    -1/sqrt(3), -1/sqrt(3),  1/sqrt(3);
     1/sqrt(3), -1/sqrt(3),  1/sqrt(3);
     1/sqrt(3),  1/sqrt(3),  1/sqrt(3);
    -1/sqrt(3),  1/sqrt(3),  1/sqrt(3)
];

gauss_weights = ones(8, 1);

% Funzioni di forma e derivate
shape_functions = @(xi, eta, zeta) [
    1/8 * (1 - xi)*(1 - eta)*(1 - zeta);
    1/8 * (1 + xi)*(1 - eta)*(1 - zeta);
    1/8 * (1 + xi)*(1 + eta)*(1 - zeta);
    1/8 * (1 - xi)*(1 + eta)*(1 - zeta);
    1/8 * (1 - xi)*(1 - eta)*(1 + zeta);
    1/8 * (1 + xi)*(1 - eta)*(1 + zeta);
    1/8 * (1 + xi)*(1 + eta)*(1 + zeta);
    1/8 * (1 - xi)*(1 + eta)*(1 + zeta)
];

shape_function_derivatives = @(xi, eta, zeta) [
    -1/8 * (1 - eta)*(1 - zeta), -1/8 * (1 - xi)*(1 - zeta), -1/8 * (1 - xi)*(1 - eta);
     1/8 * (1 - eta)*(1 - zeta), -1/8 * (1 + xi)*(1 - zeta), -1/8 * (1 + xi)*(1 - eta);
     1/8 * (1 + eta)*(1 - zeta),  1/8 * (1 + xi)*(1 - zeta), -1/8 * (1 + xi)*(1 + eta);
    -1/8 * (1 + eta)*(1 - zeta),  1/8 * (1 - xi)*(1 - zeta), -1/8 * (1 - xi)*(1 + eta);
    -1/8 * (1 - eta)*(1 + zeta), -1/8 * (1 - xi)*(1 + zeta),  1/8 * (1 - xi)*(1 - eta);
     1/8 * (1 - eta)*(1 + zeta), -1/8 * (1 + xi)*(1 + zeta),  1/8 * (1 + xi)*(1 - eta);
     1/8 * (1 + eta)*(1 + zeta),  1/8 * (1 + xi)*(1 + zeta),  1/8 * (1 + xi)*(1 + eta);
    -1/8 * (1 + eta)*(1 + zeta),  1/8 * (1 - xi)*(1 + zeta),  1/8 * (1 - xi)*(1 + eta)
];

% Inizializzazione della matrice di rigidità K
K = zeros(24, 24);

% Iterazione sui punti di Gauss
for i = 1:length(gauss_pts)
    xi = gauss_pts(i, 1);
    eta = gauss_pts(i, 2);
    zeta = gauss_pts(i, 3);
    weight = gauss_weights(i);
    
    % Calcolo delle funzioni di forma e delle loro derivate
    N = shape_functions(xi, eta, zeta);
    dN_dxi_eta_zeta = shape_function_derivatives(xi, eta, zeta);
    
    % Costruzione della matrice Jacobiana J
    J = dN_dxi_eta_zeta' * nodes;
    
    detJ = det(J);
    if detJ <= 0
        error('Jacobiano non valido in Gauss point %d', i);
    end
    
    J_inv = inv(J);
    
    % Derivata delle funzioni di forma rispetto alle coordinate globali
    dN_dx = J_inv * dN_dxi_eta_zeta';
    
    % Costruzione della matrice B
    B = zeros(6, 24);
    for j = 1:8
        B(1, 3*j-2) = dN_dx(1, j);
        B(2, 3*j-1) = dN_dx(2, j);
        B(3, 3*j)   = dN_dx(3, j);
        B(4, 3*j-2) = dN_dx(2, j);
        B(4, 3*j-1) = dN_dx(1, j);
        B(5, 3*j-2) = dN_dx(3, j);
        B(5, 3*j)   = dN_dx(1, j);
        B(6, 3*j-1) = dN_dx(3, j);
        B(6, 3*j)   = dN_dx(2, j);
    end
    
    % Calcolo della contribuzione alla matrice K
    K = K + B' * D * B * detJ * weight;
end

disp('Matrice di rigidità K:');
disp(K);