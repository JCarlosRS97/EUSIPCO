function pow_alloc = water_fill(Trans_Power,snr)
Noise_Power = 1./snr;
[S_Number,dt] = sort(Noise_Power);
for p=length(S_Number):-1:1
    T_P = (Trans_Power+sum(S_Number(1:p)))/p;
    Pt = T_P-S_Number(1:p);
    if (Pt(:)>=0)
        break
    end
end
pow_alloc = zeros(1,length(Noise_Power));
pow_alloc(dt(1:p)) = Pt;
end