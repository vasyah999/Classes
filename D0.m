function d=D0(t)
tend=0.25;
if (t<tend)
    d=0.1/tend*(tend-t);
else
    d=0;
end
end