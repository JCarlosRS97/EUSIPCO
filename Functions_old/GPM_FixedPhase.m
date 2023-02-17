function [Cout,all_steps,iter_time,myomegaunitcirc] = GPM_FixedPhase(Nt,Nr,Nris,SNR,quant_bits,Pt,Hdir,H1,H2,maxIter,Qinit,myomegaunitball,c)

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
tic
Cout = [RIScap(Hdir,H2,H1,myomegaunitcirc,Q,SNR)];
Cprev = Cout;
iter_time = toc;
while iIter<maxIter
    iIter = iIter+1;
    for iLineSearch =0:30
        %iLineSearch
        % New Q matrix
        Qnew = Q + stepsize*gracovcap(Hdir,H2,H1,myomegaunitcirc,Q,SNR);
        Qnew = cov_mat_proj_modified(Qnew,Pt*c^2);
        yunitcirc = myomegaunitcirc;
        if quant_bits > 0
             myomegaunitcircnext = quant_angle(yunitcirc,quant_bits)/c;
        else
            myomegaunitcircnext = projectontounitcircle(yunitcirc)/c;
        end
        %     objsequnitcirc(iIter) = RIScap(Hdir,H2,H1,myomegaunitcircnext,Qnew,SNR);
        Cnew = RIScap(Hdir,H2,H1,myomegaunitcircnext,Qnew,SNR);
        
%         if (((Cnew-Cprev) >= delta*(norm(Qnew-Q)^2+ ...
%                 norm(myomegaunitcircnext-myomegaunitcirc)^2)) ...
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
    all_steps(iIter) = 0;
    Cout = [Cout Cnew];
    iter_time = [iter_time toc];
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






