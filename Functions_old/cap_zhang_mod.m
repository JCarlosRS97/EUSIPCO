function [C1,Carr,alpha_opt_init,Qopt_init,iter_time] = cap_zhang_mod(Nt,Nr,Nris,SNR,L,quant_bits,Pt,Hdir,H1,H2,no_iter)

% [H1,H2,c2] = chan_mat(Nt,Nr,Nris,ht,hr,SNR,D);
Tall = {H1}; Rall = {H2}; Hall = {Hdir};
N0 = 1/SNR;
iIter = 0;
for mat_idx = 1:length(Hall)
    H = Hall{mat_idx}; T = Tall{mat_idx}; R = Rall{mat_idx};
    tic
    phases = randi([-180 180],L,Nris)*pi/180;
    alpha_all = exp(1i*phases);
    % random generated RIS phase shifts and covarinace matices 
    for k = 1:L
        alpha = alpha_all(k,:);
        Hr = H+R*diag(alpha)*T;
        Q{k} = cov_mat_opt(Nt,Nr,Hr,N0,Pt);
        C(k) = real(log2(det(eye(Nr)+(1/N0)*Hr*Q{k}*Hr')));
    end
    % choosing the best one
    [C1,kopt] = max(C);
    alpha_opt = alpha_all(kopt,:);
    Qopt = Q{kopt};
    Carr = C1; 
    % phase shift and Q matix for max. cap. 
    alpha_opt_init = alpha_opt;
    Qopt_init = Qopt;
    k = 0;
    % alternate optimization
%     while rep
   iter_time = toc;
    while iIter <no_iter
        % RIS phase optimization  
        for m = 1:Nris
            [alpha_opt_new,A{m},B{m}] = phase_shift_new(Nr,Qopt,H,T,R,Nris,m,alpha_opt,N0);
            if quant_bits>0
               alpha_opt_new(m) = quant_angle(alpha_opt_new(m),quant_bits); 
            end
            alpha_opt = alpha_opt_new;
            k = k+1;
            Hr = H+R*diag(alpha_opt)*T;
            Carr = [Carr real(log2(det(eye(Nr)+(1/N0)*Hr*Qopt*Hr')))];
            iIter = iIter+1;
            iter_time = [iter_time toc];
        end     
        Qopt = cov_mat_opt(Nt,Nr,Hr,N0,Pt);
        Hr = H+R*diag(alpha_opt_new)*T;
        Cnew = real(log2(det(eye(Nr)+(1/N0)*Hr*Qopt*Hr')));
        iIter = iIter+1;
        k = k+1;
        Carr = [Carr Cnew];
        iter_time = [iter_time toc];
    end
end
Carr = Carr(1:no_iter+1);
iter_time = iter_time(1:no_iter+1);
end