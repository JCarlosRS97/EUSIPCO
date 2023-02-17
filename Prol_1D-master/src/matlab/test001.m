
c=6.28;
matdim=200;
minEigenvalRatio = 10^-2;

[prolate_dat, iserr] = prolate_1d_crea(c,matdim, minEigenvalRatio);

figure;
semilogy((abs(prolate_dat.nu)))
ylim([10^-30,3])
title('magnitude of eigenvalues (not scaled)')

xx=linspace(-1,1,100);
prolate_ids = [1:8];
[v,dv] = prolate_1d_ev(prolate_dat, prolate_ids, xx);

figure;
plot(v(:,1:2:end));
title('even prolates')
figure;
plot(v(:,2:2:end));
title('odd prolates')

figure;
plot(dv(:,1:2:end));
title('derivative of even prolates')
figure;
plot(dv(:,2:2:end));
title('derivative of odd prolates')







