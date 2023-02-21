function [focus_f] = compute_focusing_function(RIS_pos, center, k0)

d = vecnorm(RIS_pos-center.');

focus_f = exp(-1j*k0*d).';
end