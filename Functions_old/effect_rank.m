function [N_rank] = effect_rank(eigenval)


% Computes the eigenvalues distribution
eigen_dist = eigenval./sum(abs(eigenval));

entropy_eigen = - eigen_dist.'*log(eigen_dist);

N_rank = exp(entropy_eigen);
end

