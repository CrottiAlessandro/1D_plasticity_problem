
% Define the deformation gradient F
F = [1.8, 0, 1.5;
     0, -1, 4.8;
     0, -1.6, 0];

% Calculate C = F^T * F
C = F' * F;

% Calculate eigenvalues and eigenvectors of C
[eigvecs_C, eigvals_C] = eig(C);

% Construct U
Lambda_U = diag(sqrt(diag(eigvals_C)));
U = eigvecs_C * Lambda_U * eigvecs_C';

% Calculate U^{-1}
U_inv = inv(U);

% Calculate R = F * U^{-1}
R = F * U_inv;

% Verify: R should be orthogonal
orthogonality = all(all(abs(R' * R - eye(3)) < 1e-10));

disp('Eigenvalues of C:');
disp(diag(eigvals_C));
disp('Tensor U:');
disp(U);
disp('Tensor R:');
disp(R);
disp(['R is orthogonal: ', num2str(orthogonality)]);

% Calculate B = F * F^T
B = F * F';

% Calculate eigenvalues and eigenvectors of B
[eigvecs_B, eigvals_B] = eig(B);

% Construct V
Lambda_V = diag(sqrt(diag(eigvals_B)));
V = eigvecs_B * Lambda_V * eigvecs_B';

disp('Eigenvalues of B:');
disp(diag(eigvals_B));
disp('Tensor V:');
disp(V);