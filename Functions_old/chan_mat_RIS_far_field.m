function [Hdir,H1,H2] = chan_mat_RIS_far_field(Nt,Nr,Nris,lt,lr,D,no_mat,K,f,dist_ris,varargin)

lambda = 3e8/f;     % wavelength
dt = lambda/2;      % TX antenna space
dr = lambda/2;      % RX antenna space
dris = lambda/2;    % RIS element space
k = 2*pi/lambda;    % wavenumber

if isempty(varargin)
    alpha = 1;
else
    alpha = varargin{1};
end

% geometrical placement
% x, y and z axis 
% TX array
tx_arr = zeros(3,Nt); 
tx_arr(1,:) = -dist_ris   ;
tx_arr(3,:) = lt + ((1:Nt)*dt-(Nt+1)*dt/2);
% RX array
rx_arr = zeros(3,Nr); 
rx_arr(1,:) = D-dist_ris  ;
rx_arr(3,:) =  lr + ((1:Nr)*dr-(Nr+1)*dr/2);
% RIS
RISPosition = zeros(3,Nris);
N1 = sqrt(Nris);
for i = 1:N1
    for ii = 1:N1
        RISPosition( 1, (i-1)*N1+ii ) = (i - (N1+1)/2)*dris;
        RISPosition( 2, (i-1)*N1+ii ) = (ii - (N1+1)/2)*dris;
    end
end
ris_arr = RISPosition;



% direct TX-RX paths/channel matrix
for i1 = 1:Nr
    for j1 = 1:Nt
        d(i1,j1) = norm(rx_arr(:,i1)-tx_arr(:,j1));
    end 
end
Hdir_los = exp(-1i*k*d);
tx_rx_dist = sqrt(D^2+(lt-lr)^2);
FSPL_dir = (lambda/(4*pi))^2/tx_rx_dist^alpha(1);
Hdir = Rician(Hdir_los,sqrt(FSPL_dir),no_mat,K);

% indirect paths (TX-RIS-RX)
for l1 = 1:Nris
    RISindex = rem( l1, sqrt(Nris));
    for t1 = 1:Nt
        d1(l1,t1) = dist_ris/sqrt(lt^2+dist_ris^2) * ris_arr(1,l1) + lt*(lt-tx_arr(3,t1))/sqrt(lt^2+dist_ris^2);
    end
    for r1 = 1:Nr  
%         d2(r1,l1) = -dris* sin_ref  * RISindex;
        d2(r1,l1) = -(D-dist_ris)/sqrt(lt^2+(D-dist_ris)^2) * ris_arr(1,l1) + lr*(lr-rx_arr(3,r1))/sqrt(lr^2+(D-dist_ris)^2);
    end
end



tx_ris_dist = sqrt(dist_ris^2+lt^2);
ris_rx_dist = sqrt((D-dist_ris)^2+lr^2);
Gt = 2;
Gr = 2;
G = 1;

FSPLindir = Gt*Gr* lambda^4/(256*pi^2)*...
           ((lt/tx_ris_dist*lr/ris_rx_dist))*...
           1/(tx_ris_dist*ris_rx_dist)^2;
       
%        FSPLindir = 1;
% TX-RIS
H1_los = exp(1i*k*d1);
FSPL_1 = sqrt(FSPLindir);
H1 = Rician(H1_los,FSPL_1,no_mat,K);
% RIS-RX
H2_los = exp(1i*k*d2);
FSPL_2 = 1;
H2 = Rician(H2_los,FSPL_2,no_mat,K);
end

function pos = RISPosition(N1,N2,dist,center)
d1 = (0:N1-1)-(N1-1)/2;
d2 = (0:N2-1)-(N2-1)/2;
pos{1} = center(1)+d1*dist;
pos{2} = center(2)+d2*dist;
end 

function Hout = Rician(Hlos,FSPL,no_mat,K)
Hlos = repmat(Hlos,no_mat,1);
Hnlos = sqrt(1/2)*(randn(size(Hlos))+1i*randn(size(Hlos)));
Htot = FSPL/sqrt(K+1)*(Hlos*sqrt(K)+Hnlos);
dim = size(Hlos,1)/no_mat;
for ind = 1:no_mat
   Hout{ind} = Htot((ind-1)*dim+1:ind*dim,:); 
end
end