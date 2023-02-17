function [Hdir,H1,H2] = chan_mat_RIS_new_model_ob(Nt,Nr,Nris,lt,lr,D,no_mat,K,f,dist_ris,varargin)

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
temp = [-lt,0,-dist_ris];
vector_r_t = temp/norm(temp);
tx_arr = zeros(3,Nt); 
tx_arr(1,:) = -dist_ris + ( dt*(1:Nt) -(Nt+1)*dt/2) *vector_r_t(1)  ;
tx_arr(3,:) = lt + ((1:Nt)*dt-(Nt+1)*dt/2)*vector_r_t(3);
% RX array
temp = [-lr,0,D-dist_ris];
vector_r_r = temp/norm(temp);
rx_arr = zeros(3,Nr); 
rx_arr(1,:) = D-dist_ris + ( dr*(1:Nr) -(Nr+1)*dr/2) *vector_r_r(1)  ;
rx_arr(3,:) =  lr + ((1:Nr)*dr-(Nr+1)*dr/2)*vector_r_r(3);
% RIS
RISPosition = zeros(3,Nris);
N1 = sqrt(Nris);
for i = 1:N1
    for ii = 1:N1
        RISPosition( 1, (i-1)*N1+ii ) = (i - (N1+1)/2)*dris;
        RISPosition( 2, (i-1)*N1+ii ) = (ii - (N1+1)/2)*dris;
    end
end


Constant1 = sqrt( lambda^4/(256*pi^2)) ; 

H1_los = zeros(Nris,Nt);
for i = 1:Nt
    for ii = 1:Nris
        distance_nt_nm = norm( tx_arr(:,i) - RISPosition(:,ii) );
        
        cos_theta_t_nt_nm = tx_arr(3,i)/distance_nt_nm;
        F = cos_theta_t_nt_nm;
        
        d1 = norm( tx_arr(1,i) );
        dnm = norm( RISPosition(1,ii) );
        cos_theta_tx_nt_nm = ( d1^2  + distance_nt_nm^2 - dnm^2) / 2/d1/distance_nt_nm;

        F_tx = (cos_theta_tx_nt_nm)^(Gt/2-1);
        
        Amplitude_nt_nm = Constant1*G*Gt*F_tx*F/distance_nt_nm/distance_nt_nm;
        Phase_term_nt_nm = exp( -1j*2*pi*distance_nt_nm/lambda);
        
        H1_los(ii,i) = sqrt(Amplitude_nt_nm)*Phase_term_nt_nm;
    end
end


H2_los = zeros(Nr,Nris);
for i = 1:Nr
    for ii = 1:Nris
        distance_nr_nm = norm( rx_arr(:,i) - RISPosition(:,ii) );
        
        cos_theta_r_nr_nm = rx_arr(3,i)/distance_nr_nm;
        F = cos_theta_r_nr_nm;
        
        d2 = norm( rx_arr(1,i) );
        dnm = norm( RISPosition(1,ii) );
        cos_theta_rx_nr_nm = ( d2^2  + distance_nr_nm^2 - dnm^2) / 2/d2/distance_nr_nm;
        F_rx = (cos_theta_rx_nr_nm)^(Gr/2-1);
        
        Amplitude_nr_nm = Constant1*G*Gr*F_rx*F/distance_nr_nm/distance_nr_nm;
        Phase_term_nr_nm = exp( -1j*2*pi*distance_nr_nm/lambda);
        
        H2_los(i,ii) = sqrt(Amplitude_nr_nm)*Phase_term_nr_nm;
    end
end


if isempty(varargin)
    alpha = 1;
else
    alpha = varargin{1};
end


Constant2 = sqrt( lambda*lambda/16/pi^2) ;
Hdir_los = zeros(Nr,Nt);
for i = 1:Nr
    for ii = 1:Nt
        
        vector_temp = rx_arr(:,i) - tx_arr(:,ii);
        vector1 = vector_temp/norm(vector_temp);
        vector2 = tx_arr(:,ii)/norm(tx_arr(:,ii));
        cos_theta_tx_nr_nt = - vector1'*vector2;
        if cos_theta_tx_nr_nt<0
            cos_theta_tx_nr_nt=0;
        end
        F_tx = cos_theta_tx_nr_nt^(Gt/2-1);
        
        vector3 = rx_arr(:,i)/norm(rx_arr(:,i));
        cos_theta_rx_nr_nt = vector1'*vector3;
        if cos_theta_rx_nr_nt<0
            cos_theta_rx_nr_nt=0;
        end
        F_rx = (cos_theta_rx_nr_nt)^(Gr/2-1);
        
        distance_nr_nt = norm( rx_arr(:,i) - tx_arr(:,ii) );
        
%         Amplitude_nr_nt = Constant2*Gt*Gr*F_rx*F_tx/distance_nr_nt.^4;
        Amplitude_nr_nt = (lambda/(4*pi))^2/distance_nr_nt^alpha(1);
        Phase_term_nr_nt = exp( -1j*2*pi*distance_nr_nt/lambda);
        
        Hdir_los(i,ii) = sqrt(Amplitude_nr_nt)*Phase_term_nr_nt;
    end
end


Hdir = Rician_ewise(Hdir_los,no_mat,K);
H1 = Rician_ewise(H1_los,no_mat,K);
H2 = Rician_ewise(H2_los,no_mat,K);

function Hout = Rician_ewise(Hlos,no_mat,K)
    for ind = 1:no_mat
        Hnlos = sqrt(1/2)*(randn(size(Hlos))+1i*randn(size(Hlos)));
        Htot = (Hlos*sqrt(K)+Hnlos.*Hlos)/sqrt(K+1);
        Hout{ind} = Htot;
    end
end



end
