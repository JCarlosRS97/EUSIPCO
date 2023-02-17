function [eigenfunctions,eigenvalues] = prolate_1d_gen(c,xx,prolate_ids)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

matdim = 200;
minEigenvalRatio = 1e-3;

[prolate_dat, iserr] = prolate_1d_crea(c,matdim, minEigenvalRatio);

[v,dv] = prolate_1d_ev(prolate_dat, prolate_ids, xx);

norm_factor = vecnorm(v);
v_norm = v./norm_factor;
eigenvalues_frame = prolate_dat.nu(1:length(norm_factor)).';

% Correction on the phase
theoretical_phase = 1i.^prolate_ids.';
phase_diff = angle(eigenvalues_frame) - angle(theoretical_phase);

eigenvalues = eigenvalues_frame.*exp(1j*phase_diff);
eigenfunctions = real(v_norm.*exp(-1j*phase_diff).'); % The eigenfunctions need to be real

end

