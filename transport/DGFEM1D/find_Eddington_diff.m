function [Rel_L1_diff]=find_Eddington_diff(E,Ep)

global npar dat

L1_diff=0;
E_sum=0;
for iel=1:npar.nel
    my_zone=npar.iel2zon(iel);
    %Only bother to counce differences in E if the region has a perturbed
    %sigma_t
    if dat.sigtPertRegion(my_zone)~=0 
        dx=npar.dx(iel);
        E0=E(1,iel);
        E1=E(2,iel);
        Eloc = (E1+E0)/2;
        Ep0=Ep(1,iel);
        Ep1=Ep(2,iel);
        Eploc = (Ep1+Ep0)/2;
        L1_diff=L1_diff+(abs(Eploc-Eloc)*dx);
        E_sum=E_sum+(abs(Eloc)*dx);
    end
end
Rel_L1_diff=L1_diff/E_sum;
return
    