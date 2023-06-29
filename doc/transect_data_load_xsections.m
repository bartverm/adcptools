load sample_data/rdi_muara_muntai_bend/tid.mat tid

ef(7) = EnsembleFilter;

for ce = 1:7
    ef(ce) = EnsembleFilter;
    ef(ce).bad_ensembles = tid(ce,:)==0;
end

xs = XSection(ef, mmbend);
xs(7).revert();
