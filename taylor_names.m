function par_name = taylor_names(var, dim, ord)
if ord == 0
    par_name = sprintf('%s', var);
else
    par_name = sprintf('d^%u%s/d%s^%u', ord, var, dim, ord);
end
end