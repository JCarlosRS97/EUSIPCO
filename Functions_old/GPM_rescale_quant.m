function [Cquant,quant_iter] = GPM_rescale_quant(Nt,Nr,Nris,SNR,quant_bits,Pt,Hdir,H1,H2,maxIter,Qinit,myomegaunitball,c)

% initialization
myomegaunitcirc = myomegaunitball.';
Q = Qinit; 


% the code is for minimizing -(log(det(eye(Nr)+SNR*Z*Z')) 
delta = 1e-5;
rho = 0.5;

iIter = 0;
Cprev = 0;
% res = 0;
% step_idx = 1;
stepsize = 10000; 
% start time measurement
step = 1;

tic
Cout = [RIScap(Hdir,H2,H1,myomegaunitcirc,Q,SNR)];
Cprev = Cout;
iter_time = toc;

Cquant = [Cout];
while iIter<maxIter
    iIter = iIter+1;
%      stepsize = 1000000;
    %Cprev =RIScap(Hdir,H2,H1,myomegaunitcirc,Q,SNR);
    for iLineSearch =0:30
        %iLineSearch
        % New Q matrix
        Qnew = Q + stepsize*gracovcap(Hdir,H2,H1,myomegaunitcirc,Q,SNR);
        Qnew = cov_mat_proj_modified(Qnew,Pt*c^2);
        % New RIS phase shifts
        yunitcirc = myomegaunitcirc + stepsize*gradRIScap(Hdir,H2,H1,myomegaunitcirc,Q,SNR);
        myomegaunitcircnext = projectontounitcircle(yunitcirc)/c;
        
        if (quant_bits>0) && (mod(iIter,step)==0)  
             myomegaunitcircnext_quant = quant_angles(myomegaunitcircnext,quant_bits,c);
        end
        
        Cnew = RIScap(Hdir,H2,H1,myomegaunitcircnext,Qnew,SNR);
        
        if (((Cnew-Cprev) >= delta*(norm(Qnew-Q)^2+ ...
                norm(myomegaunitcircnext-myomegaunitcirc)^2)) ...
            || (stepsize<1e-4) )
            % Update Q matrix and RIS phase shifts
            myomegaunitcirc = myomegaunitcircnext;
            Q = Qnew;
            % Update capacity
            Cprev = Cnew;
            if (mod(iIter,step)==0)               
                Cquant = [Cquant RIScap(Hdir,H2,H1,myomegaunitcircnext_quant,Qnew,SNR)];
            end
            break
        else 
            stepsize=stepsize*rho;
        end
    end
%     all_steps(iIter) = stepsize;
%     Cout = [Cout Cnew];
%     iter_time = [iter_time toc];
    quant_iter = 1:step:maxIter;
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

function y = quant_angles(x,quant_bits,c)
levels = (0:2^quant_bits-1)/(2^quant_bits);
levels = exp(1i*levels*2*pi)/c;
Nris = size(x,1);
for ind =  1:Nris
    delta = abs(x(ind,1)-levels).^2;
    [~,minind] = min(delta);
    y(ind,1) = levels(minind); 
end
end
