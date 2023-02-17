function [Hdir,H1,H2] = chan_mat_RIS_new_model_UPA_rec_ob(Nt_y,Nt_z,Nr_y,Nr_z,Nris_x,Nris_y,lt,lr,D,no_mat,Kdir,K1,K2,f,dist_ris,varargin)

lambda = 3e8/f;     % wavelength
dt = lambda/2;      % TX antenna space
dr = lambda/2;      % RX antenna space
dris = lambda/2;    % RIS element space
k = 2*pi/lambda;    % wavenumber


Gt = 2;
Gr = 2;
G = 1;



% geometrical placement
% x, y and z axis 
% TX array
Nt = Nt_y * Nt_z;
Nr = Nr_y * Nr_z;
tx_arr = zeros(3,Nt); 
% tx_arr(1,:) = -dist_ris   ;
% tx_arr(2,:) = ((1:Nt_y)*dt-(Nt_y+1)*dt/2);
% tx_arr(3,:) = lt + ((1:Nt_z)*dt-(Nt_z+1)*dt/2);
for i = 1:Nt_y
    for ii  = 1:Nt_z
        tx_arr(1,(i-1)*Nt_z+ii) = -dist_ris   ;
        tx_arr(2,(i-1)*Nt_z+ii) = ( i *dt-(Nt_y+1)*dt/2);
        tx_arr(3,(i-1)*Nt_z+ii) = lt + ( ii*dt-(Nt_z+1)*dt/2);
    end
end
% RX array
rx_arr = zeros(3,Nr); 
% rx_arr(1,:) = D-dist_ris  ;
% rx_arr(2,:) =   ((1:Nr_y)*dr-(Nr_y+1)*dr/2);
% rx_arr(3,:) =  lr + ((1:Nr)*dr-(Nr+1)*dr/2);
for i = 1:Nr_y
    for ii  = 1:Nr_z
        rx_arr(1,(i-1)*Nr_z+ii) = D-dist_ris   ;
        rx_arr(2,(i-1)*Nr_z+ii) = ( i *dr-(Nr_y+1)*dr/2);
        rx_arr(3,(i-1)*Nr_z+ii) = lr + ( ii*dr-(Nr_z+1)*dr/2);
    end
end

% RIS
Nris = Nris_x * Nris_y;
RISPosition = zeros(3,Nris);
% N1 = sqrt(Nris);
for i = 1:Nris_x
    for ii = 1:Nris_y
        RISPosition( 1, (i-1)*Nris_y+ii ) = (i - (Nris_x+1)/2)*dris;
        RISPosition( 2, (i-1)*Nris_y+ii ) = (ii - (Nris_y+1)/2)*dris;
    end
end

% Constant1 = sqrt( dris*dris*lambda*lambda/64/pi^3) ;
Constant1 = sqrt( lambda^4/(256*pi^2)) ; 
H1_los = zeros(Nris,Nt);
for i = 1:Nt
    for ii = 1:Nris
        distance_nt_nm = norm( tx_arr(:,i) - RISPosition(:,ii) );
        cos_theta_t_nt_nm = tx_arr(3,i)/distance_nt_nm;
        Amplitude_nt_nm = Constant1*G*Gt*cos_theta_t_nt_nm/distance_nt_nm/distance_nt_nm;
        Phase_term_nt_nm = exp( 1j*2*pi*distance_nt_nm/lambda);
        H1_los(ii,i) = sqrt(Amplitude_nt_nm)*Phase_term_nt_nm;
    end
end


H2_los = zeros(Nr,Nris);
for i = 1:Nr
    for ii = 1:Nris
        distance_nr_nm = norm( rx_arr(:,i) - RISPosition(:,ii) );
        cos_theta_r_nr_nm = rx_arr(3,i)/distance_nr_nm;       
        Amplitude_nr_nm = Constant1*G*Gr*cos_theta_r_nr_nm/distance_nr_nm/distance_nr_nm;
        Phase_term_nr_nm = exp( 1j*2*pi*distance_nr_nm/lambda);
        H2_los(i,ii) = sqrt(Amplitude_nr_nm)*Phase_term_nr_nm;
    end
end





% Constant2 = sqrt( lambda*lambda/16/pi^2) ;
% Hdir_los = zeros(Nr,Nt);
% for i = 1:Nr
%     for ii = 1:Nt
%         
%         vector_temp = rx_arr(:,i) - tx_arr(:,ii);
%         vector1 = vector_temp/norm(vector_temp);
%         vector2 = tx_arr(:,ii)/norm(tx_arr(:,ii));
%         cos_theta_tx_nr_nt = - vector1'*vector2;
%         if cos_theta_tx_nr_nt<0
%             cos_theta_tx_nr_nt=0;
%         end
%         F_tx = cos_theta_tx_nr_nt^(Gt/2-1);
%         
%         vector3 = rx_arr(:,i)/norm(rx_arr(:,i));
%         cos_theta_rx_nr_nt = vector1'*vector3;
%         if cos_theta_rx_nr_nt<0
%             cos_theta_rx_nr_nt=0;
%         end
%         F_rx = (cos_theta_rx_nr_nt)^(Gr/2-1);
%         
%         distance_nr_nt = norm( rx_arr(:,i) - tx_arr(:,ii) );
%         
%         Amplitude_nr_nt = Constant2*Gt*Gr*F_rx*F_tx/distance_nr_nt.^4;
% %         Amplitude_nr_nt = Constant2*Gt*Gr/distance_nr_nt/distance_nr_nt;
%         Phase_term_nr_nt = exp( -1j*2*pi*distance_nr_nt/lambda);
%         
%         Hdir_los(i,ii) = sqrt(Amplitude_nr_nt)*Phase_term_nr_nt;
%     end
% end
% if isempty(varargin)
%     alpha = 1;
% else
%     alpha = varargin{1};
% end
% direct TX-RX paths/channel matrix
for i1 = 1:Nr
    for j1 = 1:Nt
        d(i1,j1) = norm(rx_arr(:,i1)-tx_arr(:,j1));
    end 
end
% Hdir_los = exp(1i*k*d) .* sqrt( (lambda/(4*pi))^2./(d.^(3)) );
Hdir_los = exp(1i*k*d) .* sqrt( Gt*Gr*(lambda/(4*pi))^2./(d.^(3)) );
% tx_rx_dist = sqrt(D^2+(lt-lr)^2);
% FSPL_dir = (lambda/(4*pi))^2/tx_rx_dist^3;
% Hdir = Rician(Hdir_los,sqrt(FSPL_dir),no_mat,K);



% temp = max( max( abs(H2_los)));
temp =1 ;
Hdir = Rician_ewise(Hdir_los,no_mat,Kdir);
H1 = Rician_ewise(H1_los*temp,no_mat,K1);
H2 = Rician_ewise(H2_los/temp,no_mat,K2);

function Hout = Rician_ewise(Hlos,no_mat,K)
    for ind = 1:no_mat
        Hnlos = sqrt(1/2)*(randn(size(Hlos))+1i*randn(size(Hlos)));
        Htot = (Hlos*sqrt(K)+Hnlos.*Hlos)/sqrt(K+1);
        Hout{ind} = Htot;
    end
end

end
