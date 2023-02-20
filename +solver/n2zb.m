function zb = n2zb(n, m)
% Get bed elevation from n coordinate, mesh m.

n(n>max(m.nb_all)) = max(m.nb_all);
n(n<min(m.nb_all)) = min(m.nb_all);


zb = interp1(m.nb_all, m.zb_all, n, 'makima');



end