function Qnew = cov_mat_proj_modified(Qold,Pt)
    [U,D] = eig(Qold);
    Dnew = pow_proj_new(Pt,real(diag(D))).';
    Qnew = U*diag(Dnew)*U';
end 