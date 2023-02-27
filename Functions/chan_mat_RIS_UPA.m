function [Hdir,H1,H2] = chan_mat_RIS_UPA(tx_arr, rx_arr, ris_arr,f)

lambda = 3e8/f;     % wavelength
k_0 = 2*pi/lambda;    % wavenumber

%% H1 LOS
d_1 = sqrt((ris_arr(1,:).'-tx_arr(1,:)).^2 + (ris_arr(2,:).'-tx_arr(2,:)).^2 + (ris_arr(3,:).'-tx_arr(3,:)).^2);
H1_los = exp(1j*k_0*d_1)./(4*pi*d_1);



%% H2 LOS

d_2 = sqrt((rx_arr(1,:).'-ris_arr(1,:)).^2 + (rx_arr(2,:).'-ris_arr(2,:)).^2 + (rx_arr(3,:).'-ris_arr(3,:)).^2);
H2_los = exp(1j*k_0*d_2)./(4*pi*d_2);


d = sqrt((rx_arr(1,:).'-tx_arr(1,:)).^2 + (rx_arr(2,:).'-tx_arr(2,:)).^2 + (rx_arr(3,:).'-tx_arr(3,:)).^2);
Hdir_los = exp(1j*k_0*d)./(4*pi*d);

H1 = H1_los;
H2 = H2_los;
Hdir = Hdir_los;

end
