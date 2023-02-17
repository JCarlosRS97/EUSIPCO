% Script NAME:
%   
%
% DESCRIPTION:
%   
%
% INPUT:
%
% OUTPUT:
%
% ASSUMPTIONS AND LIMITATIONS:
%   
%
% REVISION HISTORY:
%   18/11/2022 - Carlos
%       * Initial implementation
%
clear, close all
addpath(genpath('Prol_1D-master/'))


% Integral norm
int2d_num = @(x,y, f) trapz(y, trapz(x, f, 2));

int_norm = @(x, y, f_1, f_2) int2d_num(x, y, f_1.*conj(f_2));


%% Configuration
lambda = 0.0857;                                                         % Wavelength
k_0 = 2*pi/lambda;
d = 20;                             % Distance to the receiving surface
N_patterns = 16;

% Transmit LIS
Delta_x_T = 0.8357;
Delta_y_T = 0.8357;
% r_T0 = [0, 0, 0]';                        % It is always centered in (0,0,0)
S_T = (2*Delta_x_T)*(2*Delta_y_T);          % Area of the transmitting surface

% Receive LIS
Delta_x_R = 0.8357;
Delta_y_R = 0.8357;
S_R = (2*Delta_x_R)*(2*Delta_y_R);          % Area of the surface
theta = 0;
phi = 0;
r_R0 = [d*sin(theta)*cos(phi),...           % Center of the receiving surface
    d*sin(theta)*sin(phi), d*cos(theta)];


% Parameters of interest
r_0 = vecnorm(r_R0);

N_sampling_x = 40;                         % Grid of the surfaces
N_sampling_y = 40;
%N_sampling_x = 3;
%N_sampling_y = 3;

%% Grid
% Generate vectors of edge
edge_Tx_vector = linspace(-Delta_x_T, Delta_x_T, N_sampling_x + 1);
edge_Ty_vector = linspace(-Delta_y_T, Delta_y_T, N_sampling_y + 1);
edge_Rx_vector = linspace(-Delta_x_R, Delta_x_R, N_sampling_x + 1);
edge_Ry_vector = linspace(-Delta_y_R, Delta_y_R, N_sampling_y + 1);

% Compute steps
step_x_R = edge_Rx_vector(2)-edge_Rx_vector(1);
step_y_R = edge_Ry_vector(2)-edge_Ry_vector(1);
step_x_T = edge_Tx_vector(2)-edge_Tx_vector(1);
step_y_T = edge_Ty_vector(2)-edge_Ty_vector(1);

% Compute the centers
x_T_vec = edge_Tx_vector(1:end-1).' + step_x_T/2;
y_T_vec = edge_Ty_vector(1:end-1).' + step_y_T/2;
[x_T_mat, y_T_mat] = meshgrid(x_T_vec, y_T_vec);

x_R_vec = edge_Rx_vector(1:end-1).' + step_x_R/2;
y_R_vec = edge_Ry_vector(1:end-1).' + step_y_R/2;
[x_R_mat, y_R_mat] = meshgrid(x_R_vec, y_R_vec);

% Position vectors
x_T = x_T_mat(:).'; y_T = y_T_mat(:).';
x_R = x_R_mat(:).'; y_R = y_R_mat(:).';

%% Direct numerical computation without approximations
in_num_pattern = zeros(N_sampling_x,N_sampling_y,N_patterns);
out_num_pattern = zeros(N_sampling_x,N_sampling_y,N_patterns);


% Compute normalization factor eq 161 tutorial Miller
delta_x_R = x_R_vec(2)-x_R_vec(1);
delta_y_R = y_R_vec(2)-y_R_vec(1);
delta_x_T = x_T_vec(2)-x_T_vec(1);
delta_y_T = y_T_vec(2)-y_T_vec(1);
norm_factor = sqrt(delta_y_T*delta_x_T*delta_y_R*delta_x_R);

% Distance between any point
r = sqrt((x_R' + r_R0(1) -x_T).^2 + (y_R'  + r_R0(2) - y_T).^2 + r_R0(3).^2);

% Compute green function
G = exp(-1j*k_0*r)./(4*pi*r);

% SVD of the sampled green function
DOF_miller = S_T*S_R/(lambda^2*r_0^2);
[U,S,V] = svds(G*norm_factor, 2*N_patterns); % To reduce computation 
eigenvalues = diag(S).^2;

% Sum of eigenvalues
sum_discrete_eigenvalues = sum(abs(G(:)*norm_factor).^2); % Sum of eigenvalues under the sampling assumption
sum_denom_approx = S_R*S_T/(4*pi*r_0)^2;        % Sum of eigenvalues under the approximation in the denominator
fprintf('Sum of eigenvalues under the sampling assumption: %f\n', sum_discrete_eigenvalues);
fprintf('Sum of eigenvalues under the approximation in the denominator: %f\n', sum_denom_approx);

for k = 1:N_patterns
    in_num_pattern(:,:,k) = reshape(V(:,k), [N_sampling_y, N_sampling_x]);
    out_num_pattern(:,:,k) = reshape(U(:,k), [N_sampling_y, N_sampling_x]);
end


%% Computing receiving function and check orthogonality

out_num_pattern = compute_out_pattern(in_num_pattern, x_T_mat, y_T_mat, x_R_mat, y_R_mat, r_R0, k_0);

SIR_in_num_db = check_ortho(x_T_vec, y_T_vec, in_num_pattern);
SIR_out_num_db = check_ortho(x_R_vec, y_R_vec, out_num_pattern);

figure(1)
stem(SIR_in_num_db), hold on, stem(SIR_out_num_db)
axis([0, 17, 0, 60])
xlabel('n')
ylabel('SIR_n')
legend('Over the transmitter', 'Over the receiver')

%% Generate analytical eigenfunctions
in_analytic_pattern = zeros(N_sampling_x,N_sampling_y,N_patterns);

% Third-order focusing function
F_T3 = @(x_T, y_T) exp(1j*k_0*((x_T.^2 + y_T.^2 - 2*r_R0(1)*x_T - 2*r_R0(2)*y_T)/(2*r_0) ...
    -(x_T.^4 + y_T.^4 - 4*r_R0(1)*x_T.^3 - 4*r_R0(2)*y_T.^3 + 2*x_T.^2.*y_T.^2 ...
    - 4*r_R0(1)*x_T.*y_T.^2 - 4*r_R0(2)*y_T.*x_T.^2 + 8*r_R0(1)*r_R0(2)*x_T.*y_T...
    + 4*r_R0(1)^2*x_T.^2 + 4*r_R0(2)^2*y_T.^2)/(8*r_0^3)));


% Prolate
omega_tx = k_0*Delta_x_R/r_0;               % Defined in the paper of Dardari and the report
omega_ty = k_0*Delta_y_R/r_0;

c_x = Delta_x_T*omega_tx*cos(theta)^2;
c_y = Delta_y_T*omega_ty;

% Conf
matdim = 1000;
minEigenvalRatio = 10^-8;
% Minimun value of the computed eigenvalues
[prolate_dat_x, iserr] = prolate_1d_crea(c_x,matdim, minEigenvalRatio);
[prolate_dat_y, iserr] = prolate_1d_crea(c_y,matdim, minEigenvalRatio);

prolate_ids = [0:10]; % Index of desired prolates
[phi_n_x,dv_x] = prolate_1d_ev(prolate_dat_x, prolate_ids, x_T_vec./Delta_x_T);
[phi_n_y,dv_y] = prolate_1d_ev(prolate_dat_y, prolate_ids, y_T_vec./Delta_y_T);
mu_n_x = abs(prolate_dat_x.nu.').^2;
mu_n_y = abs(prolate_dat_y.nu.').^2;

% Insert manually the order of the pairs (i,j)
% order_prolate = [1,1;
%     1,2;
%     2,1;
%     2,2;
%     3,1;
%     1,3;
%     3,2;
%     2,3;
%     3,3];

% % PSWF generation x
% [phi_n_x,lambda_n_x] = generatePSWF(x_T_vec./Delta_x_T,c_x);
%
% % PSWF generation y
% [phi_n_y, lambda_n_y] = generatePSWF(y_T_vec./Delta_y_T,c_y);
%
% % Computation of the self-adjoint eigenvalues
% mu_n_x = abs(lambda_n_x).^2*c_x/2/pi;
% mu_n_y = abs(lambda_n_y).^2*c_y/2/pi;

% Compute the order of importance of all the possible mode combinations
computed_modes_axis= ceil(sqrt(2*N_patterns));
combinations = combvec(mu_n_x(1:computed_modes_axis)',mu_n_y(1:computed_modes_axis)');
eigenval_prod = combinations(1,:).*combinations(2,:);
[eigenval_ord, ind] = sort(eigenval_prod, 'descend');
y_ind = ceil(ind/computed_modes_axis);
x_ind = mod(ind-1, computed_modes_axis)+1;

for k = 1 : N_patterns
    prolate_2d = phi_n_y(:,y_ind(k))*phi_n_x(:,x_ind(k))';

    %     norm_2_prol = trapz(y_T_vec,trapz(x_T_vec,prolate_2d.*conj(prolate_2d),2));
    norm_2_prol = int_norm(x_T_vec,y_T_vec, prolate_2d, prolate_2d);

    in_analytic_pattern(:,:,k) = F_T3(x_T_mat, y_T_mat).*prolate_2d/sqrt(norm_2_prol);
end

g_2 = (2*k_0*cos(theta))^-2*eigenval_ord.';


%% Computing receiving function and check orthogonality

out_analytic_pattern = compute_out_pattern(in_analytic_pattern, x_T_mat, y_T_mat, x_R_mat, y_R_mat, r_R0, k_0);

SIR_in_analytic_db = check_ortho(x_T_vec, y_T_vec, in_analytic_pattern);
SIR_out_analytic_db = check_ortho(x_R_vec, y_R_vec, out_analytic_pattern);

figure(2)
stem(SIR_in_analytic_db), hold on, stem(SIR_out_analytic_db)
axis([0, 17, 0, 60])
% title('Analytical')
xlabel('n')
ylabel('SIR_n')
legend('Over the transmitter', 'Over the receiver')

%% Pattern under third order taylor  expansion
in_3t_pattern = zeros(N_sampling_x,N_sampling_y,N_patterns);

rho = ((x_R' - x_T).^2 + (y_R' - y_T).^2 + 2*(x_R' - x_T)*r_R0(1)...
    + 2*(y_R' - y_T)*r_R0(2))/r_0^2;

r_3taylor = r_0*(1 + rho/2 - rho.^2/8);

% Compute green function
G_3taylor = exp(-1j*k_0*r_3taylor)./(4*pi*r_3taylor);

% SVD of the sampled green function
[U,S,V] = svds(G_3taylor*norm_factor, 2*N_patterns); % To reduce computation 
eigenval_3t = diag(S).^2;
for k = 1:N_patterns
    in_3t_pattern(:,:,k) = reshape(V(:,k), [N_sampling_y, N_sampling_x]);
end

%% Computing receiving function and check orthogonality

out_3t_pattern = compute_out_pattern(in_3t_pattern, x_T_mat, y_T_mat, x_R_mat, y_R_mat, r_R0, k_0);

SIR_in_3t_db = check_ortho(x_T_vec, y_T_vec, in_3t_pattern);
SIR_out_3t_db = check_ortho(x_R_vec, y_R_vec, out_3t_pattern);

figure(3)
stem(SIR_in_3t_db), hold on, stem(SIR_out_3t_db)
axis([0, 17, 0, 60])
title(' SVD 3rd order Taylor')
legend('Over the transmitter', 'Over the receiver')


%% Input pattern generation large shift
% Compute normalization factor eq 161 tutorial Miller
delta_x_R = x_R_vec(2)-x_R_vec(1);
delta_y_R = y_R_vec(2)-y_R_vec(1);
delta_x_T = x_T_vec(2)-x_T_vec(1);
delta_y_T = y_T_vec(2)-y_T_vec(1);
norm_factor = sqrt(delta_y_T*delta_x_T);

in_ls_pattern = zeros(N_sampling_x,N_sampling_y,N_patterns);

kernel = k_0^2/cos(theta)^2*sinc(omega_tx/pi*cos(theta)^2*(x_T' - x_T)).* ...
            sinc(omega_ty/pi*(y_T' - y_T));

[V,D] = eigs(norm_factor*kernel, 2*N_patterns);

for k = 1:N_patterns
    in_ls_pattern(:,:,k) =  F_T3(x_T_mat, y_T_mat).*reshape(V(:,k), [N_sampling_y, N_sampling_x]);
end

%% Computing receiving function and check orthogonality

out_ls_pattern = compute_out_pattern(in_ls_pattern, x_T_mat, y_T_mat, x_R_mat, y_R_mat, r_R0, k_0);

SIR_in_ls_db = check_ortho(x_T_vec, y_T_vec, in_ls_pattern);
SIR_out_ls_db = check_ortho(x_R_vec, y_R_vec, out_ls_pattern);

figure(4)
stem(SIR_in_ls_db), hold on, stem(SIR_out_ls_db)
axis([0, 17, 0, 60])
xlabel('n')
ylabel('SIR_n (dB)')
% title('Large shift kernel')
legend('Over the transmitter', 'Over the receiver')

%% Plot eigenvalues
DoF = S_T*S_R/r_0^2/lambda^2*cos(theta)^2;
figure()
plot(eigenvalues(1:DoF*2), 'bx'),hold on
plot(eigenval_3t(1:DoF*2), 'msquare')
plot(g_2(1:DoF*2), 'ro')
xline(DoF, 'g')
ground_truth = find(eigenvalues/eigenvalues(1)>=0.5, 1, 'last');
xline(ground_truth, 'k--'),hold off
legend('Numerical |g|^2','3rd order Taylor Approximation', 'Analytical |g|^2', 'DoF expression', 'Ground truth')
%% Plot in patterng_
N_pattern_paint = 5;
for i =1:N_pattern_paint
figure(i+10)
s = surf(x_T_mat, y_T_mat,abs(squeeze(in_num_pattern(:,:,i))));
s.EdgeColor = 'none';
view(2)
end
for i =1:N_pattern_paint
figure(i+20)
s = surf(x_T_mat, y_T_mat,abs(squeeze(in_analytic_pattern(:,:,i))));
s.EdgeColor = 'none';
view(2)
end
for i =1:N_pattern_paint
figure(i+30)
s = surf(x_T_mat, y_T_mat,abs(squeeze(in_ls_pattern(:,:,i))));
s.EdgeColor = 'none';
view(2)
end

%% Input pattern generation third order expansion
% Compute normalization factor eq 161 tutorial Miller
% delta_x_R = x_R_vec(2)-x_R_vec(1);
% delta_y_R = y_R_vec(2)-y_R_vec(1);
% delta_x_T = x_T_vec(2)-x_T_vec(1);
% delta_y_T = y_T_vec(2)-y_T_vec(1);
% norm_factor = sqrt(delta_y_T*delta_x_T);
% 
% in_k3_pattern = zeros(N_sampling_x,N_sampling_y,N_patterns);
% kernel_3t = zeros(length(x_T), length(x_T));
% length_x_T = length(x_T);
% tic
% parfor j = 1:length_x_T
%     for i = 1:length_x_T
%         x_T1 = x_T(i); y_T1 = y_T(i);
%         x_T2 = x_T(j); y_T2 = y_T(j);
%         D_integrand_3 = exp(-1j*k_0*((x_R_mat.*(x_T1 - x_T2) + y_R_mat.*(y_T1 - y_T2))/r_0...
%             + (-2*x_R_mat.*(x_T1.^3 - x_T2.^3) - 2*y_R_mat.*(y_T1.^3 - y_T2.^3)...
%             - 2*x_R_mat.*(x_T1.*y_T1.^2 - x_T2.*y_T2.^2) - 2*y_R_mat.*(x_T1.^2.*y_T1 - x_T2.^2.*y_T2)...
%             + 4*(r_R0(1)*y_R_mat + x_R_mat.*y_R_mat + r_R0(2)*x_R_mat)*(x_T1.*y_T1 - x_T2.*y_T2)...
%             + (6*r_R0(1)*x_R_mat + 3*x_R_mat.^2 + y_R_mat.^2 + 2*r_R0(2)*y_R_mat)*(x_T1.^2 - x_T2.^2) ...
%             + (x_R_mat.^2 + 2*r_R0(1)*x_R_mat + 6*r_R0(2)*y_R_mat + 3*y_R_mat.^2)*(y_T1.^2 - y_T2.^2)...
%             - 2*(2*r_R0(1)^2*x_R_mat + 3*r_R0(1)*x_R_mat.^2 + r_R0(1)*y_R_mat.^2 + 2*r_R0(2)*r_R0(1)*y_R_mat + x_R_mat.^3 ...
%             + x_R_mat.*y_R_mat.^2 + 2*r_R0(2)*x_R_mat.*y_R_mat).*(x_T1 - x_T2)...
%             -2*(x_R_mat.^2*r_R0(2) + x_R_mat.^2.*y_R_mat + 2*r_R0(1)*x_R_mat*r_R0(2) + 2*r_R0(1)*x_R_mat.*y_R_mat + 2*r_R0(2)^2*y_R_mat ...
%             + 3*r_R0(2)*y_R_mat.^2 + y_R_mat.^3).*(y_T1 - y_T2))/(4*r_0^3)));
% 
%         kernel_3t(j, i) = int2d_num(x_R_vec, y_R_vec, D_integrand_3);
%     end
% end
% toc
% 
% [V,D] = eigs(norm_factor*kernel_3t, N_patterns);
% 
% for k = 1:N_patterns
%     in_k3_pattern(:,:,k) =  F_T3(x_T_mat, y_T_mat).*reshape(V(:,k), [N_sampling_y, N_sampling_x]);
% end

% %% Computing receiving function and check orthogonality
% 
% out_k3_pattern = compute_out_pattern(in_k3_pattern, x_T_mat, y_T_mat, x_R_mat, y_R_mat, r_R0, k_0);
% 
% SIR_in_k3_db = check_ortho(x_T_vec, y_T_vec, in_k3_pattern);
% SIR_out_k3_db = check_ortho(x_R_vec, y_R_vec, out_k3_pattern);
% 
% figure(5)
% stem(SIR_in_k3_db), hold on, stem(SIR_out_k3_db)
% axis([0, 17, 0, 60])
% title('Kernel 3rd order')

%% Save fig
% exportgraphics(figure(2),'Eigenfunctions\SIR_analytic_oneshift.pdf','ContentType','vector')
% exportgraphics(figure(1),'Eigenfunctions\SIR_num_oneshift.pdf','ContentType','vector')
% exportgraphics(figure(4),'Eigenfunctions\SIR_int_oneshift.pdf','ContentType','vector')
% 
% exportgraphics(figure(11),'Eigenfunctions\num_mode1_oneshift.pdf','ContentType','image')
% exportgraphics(figure(21),'Eigenfunctions\Analytic_mode1_oneshift.pdf','ContentType','image')
% exportgraphics(figure(31),'Eigenfunctions\Int_mode1_oneshift.pdf','ContentType','image')
% exportgraphics(figure(11),'Eigenfunctions\num_mode2_oneshift.pdf','ContentType','image')
% exportgraphics(figure(21),'Eigenfunctions\Analytic_mode2_oneshift.pdf','ContentType','image')
% exportgraphics(figure(31),'Eigenfunctions\Int_mode2_oneshift.pdf','ContentType','image')
% exportgraphics(figure(12),'Eigenfunctions\num_mode2_oneshift.pdf','ContentType','image')
% exportgraphics(figure(22),'Eigenfunctions\Analytic_mode2_oneshift.pdf','ContentType','image')
% exportgraphics(figure(32),'Eigenfunctions\Int_mode2_oneshift.pdf','ContentType','image')
% exportgraphics(figure(13),'Eigenfunctions\num_mode3_oneshift.pdf','ContentType','image')
% exportgraphics(figure(23),'Eigenfunctions\Analytic_mode3_oneshift.pdf','ContentType','image')
% exportgraphics(figure(33),'Eigenfunctions\Int_mode3_oneshift.pdf','ContentType','image')
% exportgraphics(figure(14),'Eigenfunctions\num_mode4_oneshift.pdf','ContentType','image')
% exportgraphics(figure(24),'Eigenfunctions\Analytic_mode4_oneshift.pdf','ContentType','image')
% exportgraphics(figure(34),'Eigenfunctions\Int_mode4_oneshift.pdf','ContentType','image')
% exportgraphics(figure(15),'Eigenfunctions\num_mode5_oneshift.pdf','ContentType','image')
% exportgraphics(figure(25),'Eigenfunctions\Analytic_mode5_oneshift.pdf','ContentType','image')
% exportgraphics(figure(35),'Eigenfunctions\Int_mode5_oneshift.pdf','ContentType','image')
%% Functions

function out_pattern = compute_out_pattern(in_pattern, x_T_mat, y_T_mat, x_R_mat, y_R_mat, r_R0, k_0)

% Integral norm
int2d_num = @(x,y, f) trapz(y, trapz(x, f, 2));
int_norm = @(x, y, f_1, f_2) int2d_num(x, y, f_1.*conj(f_2));

N_sampling_x = size(in_pattern,1);
N_sampling_y = size(in_pattern,2);
N_patterns = size(in_pattern,3);

x_T_vec = x_T_mat(1,:);
y_T_vec = y_T_mat(:,1).';


out_pattern = zeros(N_sampling_x,N_sampling_y,N_patterns);

parfor k = 1:N_patterns
    k_pattern = zeros(N_sampling_y, N_sampling_x);
    for i = 1:N_sampling_x
        for j = 1:N_sampling_y
            % Distance between any point
            r = sqrt((x_R_mat(j,i) + r_R0(1) -x_T_mat).^2 + (y_R_mat(j,i)  + r_R0(2) - y_T_mat).^2 + r_R0(3).^2);

            % Compute green function
            G = exp(-1j*k_0*r)./(4*pi*r);

            integrand = G.*squeeze(in_pattern(:,:,k));
            k_pattern(j,i) = int2d_num(x_T_vec, y_T_vec, integrand);
        end
    end
    out_pattern(:,:,k) = k_pattern;
end

end

function [SIR_db] = check_ortho(x_vec, y_vec, pattern)
% Integral norm
int2d_num = @(x,y, f) trapz(y, trapz(x, f, 2));
int_norm = @(x, y, f_1, f_2) int2d_num(x, y, f_1.*conj(f_2));

N_patterns = size(pattern,3);

ortho_matrix = zeros(N_patterns,N_patterns);

for i = 1:N_patterns
    for j = 1:N_patterns
        ortho_matrix(i,j) = int_norm(x_vec,y_vec, squeeze(pattern(:,:,i)), squeeze(pattern(:,:,j)));
    end
end

inter_matrix = abs(ortho_matrix).^2;

SIR = diag(inter_matrix)./(sum(inter_matrix,2)-diag(inter_matrix));

SIR_db = 10*log10(SIR);
end