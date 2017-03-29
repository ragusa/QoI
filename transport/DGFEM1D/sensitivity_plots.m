function sensitivity_plots
x = -0.1:0.005:0.1;
y = -0.1:0.005:0.1;

[xmesh,ymesh]=meshgrid(x,y);
z=arrayfun(@findError,xmesh,ymesh);
surf(xmesh,ymesh,z)
end

function val=findError(sigsFactor,sigtFactor)
    [pqoi_sn_f,pqoi_vef_math_adj] = sn1d_iterator(sigsFactor,sigtFactor);
    val=(pqoi_sn_f - pqoi_vef_math_adj)/pqoi_sn_f;
end














