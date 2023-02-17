function Q = cov_mat_opt(Nt,Nr,H,N0,Pt)

D = min(Nt,Nr);
[~,Lambda,V] = svd(H);
Lambda = Lambda(1:D,1:D);
V = V(:,1:D);
snr = diag(Lambda).^2/N0;
pow_alloc = water_fill(Pt,snr);
Q = V*diag(pow_alloc)*V';
end 