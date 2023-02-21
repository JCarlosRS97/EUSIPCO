function [ris_arr1,ris_arr2]  = place_antennas(N_1, N_2,  d_e)
ris_v_1 = ((1:N_1)*d_e-(N_1+1)*d_e/2).';
ris_v_2 = ((1:N_2)*d_e-(N_2+1)*d_e/2).';

[ris_arr_1, ris_arr_2] = meshgrid(ris_v_1, ris_v_2);

ris_arr1 = ris_arr_1(:);
ris_arr2 = ris_arr_2(:);

end