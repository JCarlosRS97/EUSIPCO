function [alpha_out,A,B] = phase_shift_new(Nr,Q,H,T,R,M,m,alpha,N0)

[Uq,Sigmaq] = eig(Q);
Hprim = H*Uq*(Sigmaq^(1/2));
Tprim = T*Uq*(Sigmaq^(1/2));
Mind = 1:M;
Mrest = Mind(Mind~=m);
% T1 = T(Mrest,:); R1 = R(:,Mrest);
% alpha1 = alpha(Mrest);
S = Hprim;
Tprim1 = Tprim'; 
rm = R(:,m); tm = Tprim1(:,m);
for k = 1:M-1
    ind = Mrest(k);
    S = S+alpha(ind)*R(:,ind)*Tprim1(:,ind)';
end
A = eye(Nr)+(1/N0)*(S*S')+(1/N0)*rm*(tm')*tm*(rm');
B = (1/N0)*rm*(tm')*S';
alpha_out = alpha;
if trace(inv(A)*B)~=0
    Sigmam = eig(inv(A)*B).';
    [~,idx] = max(abs(Sigmam));
    alpha_out(m) = exp(-1i*phase(Sigmam(idx)));
else
    alpha_out(m) = 1;
end
end