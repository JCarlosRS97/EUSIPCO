function [Cout,Q, myomegaunitcirc] = GPM_FixedPhase_Q(Nt,Nr,Nris,SNR,quant_bits,Pt,Hdir,H1,H2,maxIter,Qinit,myomegaunitball,c)

% random phase shift initialization for RIS elements 
% myomegaunitball = exp(-1i*rand(Nris,1)*2*pi);
% myomegaunitball = ones(1,Nris)'; % exp(-1i*rand(Nris,1)*2*pi);
myomegaunitcirc = myomegaunitball.';
% Qinit = eye(Nt)*(Pt/Nt);
Q = Qinit; 


rho = 0.5;

iIter = 0;
Cprev = 0;


delta = 1e-5;
stepsize = 1e7; 
% start time measurement
Cout = [RIScap(Hdir,H2,H1,myomegaunitcirc,Q,SNR)];
Cprev = Cout;
while iIter<maxIter
    iIter = iIter+1;
    Qgrad = gracovcap(Hdir,H2,H1,myomegaunitcirc,Q,SNR);
    for iLineSearch =0:30
        %iLineSearch
        % New Q matrix
        Qnew = Q + stepsize*Qgrad;
        Qnew = cov_mat_proj_modified(Qnew,Pt*c^2);

        %     objsequnitcirc(iIter) = RIScap(Hdir,H2,H1,myomegaunitcircnext,Qnew,SNR);
        Cnew = RIScap(Hdir,H2,H1,myomegaunitcirc,Qnew,SNR);
        
%         if (((Cnew-Cprev) >= delta*(norm(Qnew-Q)^2+ ...
%                 norm(myomegaunitcircnext-myomegaunitcirc)^2)) ...
        if (((Cnew-Cprev) >= delta*norm(Qnew-Q)^2) ...
            || (stepsize<1e-20) )
            % Update Q matrix
            Q = Qnew;
            % Update capacity
            Cprev = Cnew;
            break
        else 
            stepsize=stepsize*rho;
        end
    end
    Cout = [Cout Cprev];
    if (iIter > 30 && (Cout(iIter)/Cout(iIter-20)-1 <1e-4))
        break
    end
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






