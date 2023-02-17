function [Cout,all_steps1,all_steps2,iter_time,myomegaunitcirc] = GPM_Stefan_carlos(Nt,Nr,Nris,SNR,quant_bits,Pt,Hdir,H1,H2,maxIter,Qinit,myomegaunitball,c)

% random phase shift initialization for RIS elements 
% myomegaunitball = exp(-1i*rand(Nris,1)*2*pi);
% myomegaunitball = ones(1,Nris)'; % exp(-1i*rand(Nris,1)*2*pi);
myomegaunitcirc = myomegaunitball.';
% Qinit = eye(Nt)*(Pt/Nt);
Q = Qinit; 

% % un-optimized capacity 
% Z = sqrt(SNR)*(Hdir+H2*diag(myomegaunitcirc)*H1);
% C = real(log(det(eye(Nr)+Z*Q*Z')))/log(2);
% C = RIScap(Hdir,H2,H1,myomegaunitcirc,Q,SNR);

% stepsize_val = [100 70 50 40 30 20 10 5 2 1 0.1 1e-2 1e-3 1e-4];

% maxIter = 1000;
% objsequnitcirc = zeros(1,maxIter);


% the code is for minimizing -(log(det(eye(Nr)+SNR*Z*Z')) 
% delta = 1e-5;

rho1 = 0.5;
rho2 = 0.5;

iIter = 0;
Cprev = 0;
% res = 0;
% step_idx = 1;
% stepsize = 10000; 

delta1 = 1e-1/length(myomegaunitcirc);
stepsize1 = 1e7; 
delta2 = 1e-3/length(Q(:));
stepsize2 = 1e5;
% start time measurement
tic
Cout = [RIScap(Hdir,H2,H1,myomegaunitcirc,Q,SNR)];
Cprev = Cout;
iter_time = toc;
while iIter<maxIter
    iIter = iIter+1;
%      stepsize = 1000000;
    %Cprev =RIScap(Hdir,H2,H1,myomegaunitcirc,Q,SNR);

    % Opt RIS
    PhiGrad = gradRIScap(Hdir,H2,H1,myomegaunitcirc,Q,SNR);
    for iLineSearch =0:30
        %iLineSearch
        % New RIS phase shifts
        yunitcirc = myomegaunitcirc + stepsize1*PhiGrad;
        %norm(gradRIScap(Hdir,H2,H1,myomegaunitcirc,Q,SNR))
        if quant_bits > 0
             myomegaunitcircnext = quant_angle(yunitcirc,quant_bits)/c;
        else
            myomegaunitcircnext = projectontounitcircle(yunitcirc)/c;
        end
        %     objsequnitcirc(iIter) = RIScap(Hdir,H2,H1,myomegaunitcircnext,Qnew,SNR);
        Cnew = RIScap(Hdir,H2,H1,myomegaunitcircnext,Q,SNR);
        
%         if (((Cnew-Cprev) >= delta*(norm(Qnew-Q)^2+ ...
%                 norm(myomegaunitcircnext-myomegaunitcirc)^2)) ...
        if (((Cnew-Cprev) >= delta1* ...
                norm(myomegaunitcircnext-myomegaunitcirc)^2) ...
            || (stepsize1<1e-4) )
            % Update Q matrix and RIS phase shifts
            myomegaunitcirc = myomegaunitcircnext;
%             Q = Qnew;
            % Update capacity
            Cprev = Cnew;
            break
        else 
            stepsize1=stepsize1*rho1;
        end
    end

    % Opt Q
    Qgrad = gracovcap(Hdir,H2,H1,myomegaunitcirc,Q,SNR);
    for iLineSearch =0:30
        %iLineSearch
        % New Q matrix
        Qnew = Q + stepsize2*Qgrad;
        %norm(gracovcap(Hdir,H2,H1,myomegaunitcirc,Q,SNR))
        Qnew = cov_mat_proj_modified(Qnew,Pt*c^2);

        %     objsequnitcirc(iIter) = RIScap(Hdir,H2,H1,myomegaunitcircnext,Qnew,SNR);
        Cnew = RIScap(Hdir,H2,H1,myomegaunitcircnext,Qnew,SNR);
        
%         if (((Cnew-Cprev) >= delta*(norm(Qnew-Q)^2+ ...
%                 norm(myomegaunitcircnext-myomegaunitcirc)^2)) ...
        if (((Cnew-Cprev) >= delta2*norm(Qnew-Q)^2) ...
            || (stepsize2<1e-20) )
            % Update Q matrix and RIS phase shifts
            Q = Qnew;
            % Update capacity
            Cprev = Cnew;
            break
        else 
            stepsize2=stepsize2*rho2;
        end
    end

    all_steps1(iIter) = stepsize1;
    all_steps2(iIter) = stepsize2;
    Cout = [Cout Cprev];
    iter_time = [iter_time toc];

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






