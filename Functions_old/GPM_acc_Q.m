function Cout = GPM_acc_Q(Nt,Nr,Nris,ht,hr,SNR,D,quant_bits,Pt,Hdir,H1,H2,maxIter,Qinit,myomegaunitball,c)

% random phase shift initialization for RIS elements 
% myomegaunitball = exp(-1i*rand(Nris,1)*2*pi);
% myomegaunitball = ones(1,Nris)'; % exp(-1i*rand(Nris,1)*2*pi);
myomegaunitball = myomegaunitball.';
myomegaunitcirc = myomegaunitball;
% Qinit = eye(Nt)*(Pt/Nt);

Q = Qinit; 

% un-optimized capacity 
Z = sqrt(SNR)*(Hdir+H2*diag(myomegaunitcirc)*H1);
C = real(log(det(eye(Nr)+Z*Q*Z')))/log(2);
C = RIScap(Hdir,H2,H1,myomegaunitcirc,Q,SNR);


% stepsize_val = [70 50 40 30 20 10 5 2 1 0.1 0.01 0.001];
% stepsize_val = stepsize_val/sqrt(SNR);
% maxIter = 1000;
% objsequnitcirc = zeros(1,maxIter);
stepsize1 = 10000;
stepsize2 = 10000;
delta = 1e-5;
rho = 0.5;
% the code is for minimizing -(log(det(eye(Nr)+SNR*Z*Z')) 
iIter = 0;

ThetaN = myomegaunitball;
ThetaN_prev = ThetaN;
ThetaN_bar = ThetaN;
QN = Qinit;
QN_prev = QN;
QN_bar = QN;
tn = 1; tn_prev = 1;

C_bar_prev = 0; 
C_dot_prev = 0; 
Cout = [RIScap(Hdir,H2,H1,myomegaunitcirc,Q,SNR)];
while iIter<maxIter
    
    iIter = iIter+1;
    % Steps 3 and 4 
%     ThetaN_tilde = ThetaN + tn_prev/tn*(ThetaN_bar-ThetaN) + ...
%         (tn_prev-1)/tn*(ThetaN-ThetaN_prev);
    QN_tilde = QN + tn_prev/tn*(QN_bar-QN) + ...
        (tn_prev-1)/tn*(QN-QN_prev);
    C_bar_prev = RIScap(Hdir,H2,H1,ThetaN,QN_tilde,SNR);
%     stepsize1 = 1e6;
    for iLineSearch =0:30 % for Steps 5 and 6
%         ThetaN_bar = ThetaN_tilde + stepsize1*gradRIScap(Hdir,H2,H1,ThetaN_tilde,QN_tilde,SNR);
%         ThetaN_bar = projectontounitcircle(ThetaN_bar)/c;
        QN_bar = QN_tilde + stepsize1*gracovcap(Hdir,H2,H1,ThetaN,QN_tilde,SNR);
        QN_bar = cov_mat_proj_modified(QN_bar,Pt*c^2);
        C_bar_new = RIScap(Hdir,H2,H1,ThetaN_bar,QN_bar,SNR);
        if((C_bar_new-C_bar_prev)< delta*(norm(QN_bar-QN_tilde)^2))
           stepsize1=stepsize1*rho;
           if(stepsize1<1e-4)
               break
           end
        else
            break; % break the line search procedure
        end
    end
    C_dot_prev = RIScap(Hdir,H2,H1,ThetaN,QN,SNR);
%     stepsize2 = 1e6;
    for iLineSearch=0:30
            % Steps 7 and 8
%         ThetaN_dot = ThetaN + stepsize2*gradRIScap(Hdir,H2,H1,ThetaN,QN,SNR);
%         ThetaN_dot = projectontounitcircle(ThetaN_dot)/c;
        QN_dot = QN + stepsize2*gracovcap(Hdir,H2,H1,ThetaN,QN,SNR);
        QN_dot = cov_mat_proj_modified(QN_dot,Pt*c^2);
        C_dot_new = RIScap(Hdir,H2,H1,ThetaN,QN_dot,SNR);
        if((C_dot_new-C_dot_prev) < delta*(norm(QN_dot-QN)^2))
           stepsize2=stepsize2*rho;
           if(stepsize2<1e-4)
               break
           end
        else
            break; % break the line search procedure
        end
    end
        
    % Step 9
    
    if C_bar_new>C_dot_new
%         ThetaN = ThetaN_bar;
        QN = QN_bar;
    else
%         ThetaN = ThetaN_dot;
        QN = QN_dot;
    end


    % Step 10
    tn_prev = tn;
    tn = (sqrt(4*tn^2+1)+1)/2;

    % Update capacity
    Cout = [Cout max(C_bar_new,C_dot_new)];
    
end
end



function y = gracovcap(Hdir,H2,H1,myomega,Q,SNR)
Z = sqrt(SNR)*(Hdir+H2*diag(myomega)*H1);
Nr = size(H2,1);
y = Z'*inv(eye(Nr)+Z*Q*Z')*Z;
end

function y = gradRIScap(Hdir,H2,H1,myomega,Q,SNR)
Z = sqrt(SNR)*(Hdir+H2*diag(myomega)*H1);
Nr = size(H2,1);
y = sqrt(SNR)*diag(H2'*inv(eye(Nr)+Z*Q*Z')*Z*Q*H1');
end

function y = RIScap(Hdir,H2,H1,myomega,Q,SNR)
Z = sqrt(SNR)*(Hdir+H2*diag(myomega)*H1);
Nr = size(H2,1);
y = real(log(det(eye(Nr)+Z*Q*Z')))/log(2);
end
