function [pswf, eigenval_ord] = compute_analytical_pattern(c_x, c_y, x,y, k_0, N_patterns)
% conf
matdim = 1000;
minEigenvalRatio = 10^-40;

% Minimun value of the computed eigenvalues
[prolate_dat_x, iserr] = prolate_1d_crea(c_x,matdim, minEigenvalRatio);
[prolate_dat_y, iserr] = prolate_1d_crea(c_y,matdim, minEigenvalRatio);

prolate_ids = [0:10]; % Index of desired prolates
[phi_n_x,dv_x] = prolate_1d_ev(prolate_dat_x, prolate_ids, x);
[phi_n_y,dv_y] = prolate_1d_ev(prolate_dat_y, prolate_ids, y);
mu_n_x = abs(prolate_dat_x.nu.').^2;
mu_n_y = abs(prolate_dat_y.nu.').^2;


% Compute the order of importance of all the possible mode combinations
computed_modes_axis= ceil(sqrt(2*N_patterns));
combinations = combvec(mu_n_x(1:computed_modes_axis)',mu_n_y(1:computed_modes_axis)');
eigenval_prod = combinations(1,:).*combinations(2,:);
[eigenval_ord, ind] = sort(eigenval_prod, 'descend');
y_ind = ceil(ind/computed_modes_axis);
x_ind = mod(ind-1, computed_modes_axis)+1;

pswf = zeros(length(x)*length(y));
for k = 1 : N_patterns
    prolate_2d = phi_n_y(:,y_ind(k))*phi_n_x(:,x_ind(k))';
    norm_2_prol = sqrt(sum(abs(prolate_2d(:)).^2));

    pswf(:,k) = prolate_2d(:)/norm_2_prol;
end


end