function [Cout] = GPM_rescale_CSI(Nt,Nr,Nris,SNR,Pt,Hdir,H1,H2,Hdir_err,H1_err,H2_err,maxIter,Qinit,myomegaunitball,c)

myomegaunitcirc = myomegaunitball.';
Q = Qinit; 

delta = 1e-4;
rho = 0.5;

iIter = 0;
stepsize = 10000; 
% start time measurement
Cout = [RIScap(Hdir,H2,H1,myomegaunitcirc,Q,SNR)];
Cprev = Cout;

while iIter<maxIter
    iIter = iIter+1;
    for iLineSearch =0:30
        %iLineSearch
        % New Q matrix
        Qnew = Q + stepsize*gracovcap(Hdir_err,H2_err,H1_err,myomegaunitcirc,Q,SNR);
        Qnew = cov_mat_proj_modified(Qnew,Pt*c^2);
        % New RIS phase shifts
        yunitcirc = myomegaunitcirc + stepsize*gradRIScap(Hdir_err,H2_err,H1_err,myomegaunitcirc,Q,SNR);
        myomegaunitcircnext = projectontounitcircle(yunitcirc)/c;
 
        Cnew = RIScap(Hdir,H2,H1,myomegaunitcircnext,Qnew,SNR);
      
        if (((Cnew-Cprev) >= delta*(norm(Qnew-Q)^2+ ...
                norm(myomegaunitcircnext-myomegaunitcirc)^2)) ...
            || (stepsize<1e-4) )
            % Update Q matrix and RIS phase shifts
            myomegaunitcirc = myomegaunitcircnext;
            Q = Qnew;
            % Update capacity
            Cprev = Cnew;
            break
        else 
            stepsize=stepsize*rho;
        end
    end
    Cout = [Cout Cnew];
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



