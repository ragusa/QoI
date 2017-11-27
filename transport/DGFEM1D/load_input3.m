function [nel_zone,width_zone,sigt,sigs,qvf,incf,qva,inca]=load_input3(pb_ID)
% load XS (total and isotropic scattering) for each material
% and load each volumetric external (isotropic) source intensity
% non-isotropic data is not considered here

global dat snq npar

sn=snq.n_dir;
mu=snq.mu;

%%%%%%%%%%%%%%%%%%%%%% select input data

switch pb_ID
case 84 %
        % number of elements per zone
        nel_zone = [ 10 10 10 10 10 10]*4;
        % width of each zone
        width_zone = [1 1 1 1 1 1];
        % sigt/sigs per zone
        sigt=[1 1 1 1e-8 1e-8 1e-8];
        sigs=[1 1 1 0 0 0];
        % volumetric source value, per zone
        qvf=[1 1 1 0 0 0];
        % incoming flux values
        incf(1:sn) = 0;
        incf((sn/2)+1:sn) = 1;
        % volumetric response value, per zone
        qva=[0 0 1 1 0 0];
        % incoming adj flux values
        inca(1:sn) = 0;
        %Regions to be perturbed. Use value of 1 to specify
        dat.sigaPertRegion=[0 0 0 0 0 0];
        dat.sigsPertRegion=[1 1 1 0 0 0];
        dat.sourcePertRegion=[1 1 1 1 1 1];
        dat.incPertRegion(1:sn) = 1;
        
case 85 %
        % number of elements per zone
        nel_zone = [ 10 10 10 10 10 10]*4;
        % width of each zone
        width_zone = [1 1 1 1 1 1];
        % sigt/sigs per zone
        sigt=[1 1 1 1e-8 1e-8 1e-8];
        sigs=[1 1 1 0 0 0];
        % volumetric source value, per zone
        qvf=[1 1 1 0 0 0];
        % incoming flux values
        incf(1:sn) = 0;
        incf((sn/2)+1:sn) = 1;
        % volumetric response value, per zone
        qva=[0 0 1 0 0 0];
        % incoming adj flux values
        inca(1:sn) = 0;
        %Regions to be perturbed. Use value of 1 to specify
        dat.sigaPertRegion=[0 0 0 0 0 0];
        dat.sigsPertRegion=[1 1 1 0 0 0];
        dat.sourcePertRegion=[1 1 1 1 1 1];
        dat.incPertRegion(1:sn) = 1;
       
case 86 %
        % number of elements per zone
        nel_zone = [ 10 10 10 10 10 10]*4;
        % width of each zone
        width_zone = [1 1 1 1 1 1];
        % sigt/sigs per zone
        sigt=[1 1 1 1e-8 1e-8 1e-8];
        sigs=[1 1 1 0 0 0];
        % volumetric source value, per zone
        qvf=[1 1 1 0 0 0];
        % incoming flux values
        incf(1:sn) = 0;
        incf((sn/2)+1:sn) = 1;
        % volumetric response value, per zone
        qva=[0 0 0 1 0 0];
        % incoming adj flux values
        inca(1:sn) = 0;
        %Regions to be perturbed. Use value of 1 to specify
        dat.sigaPertRegion=[0 0 0 0 0 0];
        dat.sigsPertRegion=[1 1 1 0 0 0];
        dat.sourcePertRegion=[1 1 1 1 1 1];
        dat.incPertRegion(1:sn) = 1;   
        
case 87 %
        % number of elements per zone
        nel_zone = [ 10 10 10 10 10 10]*4;
        % width of each zone
        width_zone = [1 1 1 1 1 1];
        % sigt/sigs per zone
        sigt=[1 1 1 1 1 1];
        sigs=[1 1 1 0 0 0];
        % volumetric source value, per zone
        qvf=[1 1 1 0 0 0];
        % incoming flux values
        incf(1:sn) = 0;
        incf((sn/2)+1:sn) = 1;
        % volumetric response value, per zone
        qva=[0 0 1 1 0 0];
        % incoming adj flux values
        inca(1:sn) = 0;
        %Regions to be perturbed. Use value of 1 to specify
        dat.sigaPertRegion=[0 0 0 1 1 1];
        dat.sigsPertRegion=[1 1 1 0 0 0];
        dat.sourcePertRegion=[1 1 1 1 1 1];
        dat.incPertRegion(1:sn) = 1; 
        
case 88 %
        % number of elements per zone
        nel_zone = [ 10 10 10 10 10 10]*4;
        % width of each zone
        width_zone = [1 1 1 1 1 1];
        % sigt/sigs per zone
        sigt=[1 1 1 1 1 1];
        sigs=[1 1 1 0 0 0];
        % volumetric source value, per zone
        qvf=[1 1 1 0 0 0];
        % incoming flux values
        incf(1:sn) = 0;
        incf((sn/2)+1:sn) = 1;
        % volumetric response value, per zone
        qva=[0 0 1 0 0 0];
        % incoming adj flux values
        inca(1:sn) = 0;
        %Regions to be perturbed. Use value of 1 to specify
        dat.sigaPertRegion=[0 0 0 1 1 1];
        dat.sigsPertRegion=[1 1 1 0 0 0];
        dat.sourcePertRegion=[1 1 1 1 1 1];
        dat.incPertRegion(1:sn) = 1; 
        
case 89 %
        % number of elements per zone
        nel_zone = [ 10 10 10 10 10 10]*4;
        % width of each zone
        width_zone = [1 1 1 1 1 1];
        % sigt/sigs per zone
        sigt=[1 1 1 1 1 1];
        sigs=[1 1 1 0 0 0];
        % volumetric source value, per zone
        qvf=[1 1 1 0 0 0];
        % incoming flux values
        incf(1:sn) = 0;
        incf((sn/2)+1:sn) = 1;
        % volumetric response value, per zone
        qva=[0 0 0 1 0 0];
        % incoming adj flux values
        inca(1:sn) = 0;
        %Regions to be perturbed. Use value of 1 to specify
        dat.sigaPertRegion=[0 0 0 1 1 1];
        dat.sigsPertRegion=[1 1 1 0 0 0];
        dat.sourcePertRegion=[1 1 1 1 1 1];
        dat.incPertRegion(1:sn) = 1;
        
case 90 %
        % number of elements per zone
        nel_zone = [ 10 10 10 10 10 10]*4;
        % width of each zone
        width_zone = [1 1 1 1 1 1];
        % sigt/sigs per zone
        sigt=[1 1 1 1 1 1];
        sigs=[1 1 1 1 1 1];
        % volumetric source value, per zone
        qvf=[1 1 1 0 0 0];
        % incoming flux values
        incf(1:sn) = 0;
        incf((sn/2)+1:sn) = 1;
        % volumetric response value, per zone
        qva=[0 0 1 1 0 0];
        % incoming adj flux values
        inca(1:sn) = 0;
        %Regions to be perturbed. Use value of 1 to specify
        dat.sigaPertRegion=[0 0 0 0 0 0];
        dat.sigsPertRegion=[1 1 1 0 0 0];
        dat.sourcePertRegion=[1 1 1 1 1 1];
        dat.incPertRegion(1:sn) = 1;
        
case 91 %
        % number of elements per zone
        nel_zone = [ 10 10 10 10 10 10]*4;
        % width of each zone
        width_zone = [1 1 1 1 1 1];
        % sigt/sigs per zone
        sigt=[1 1 1 1 1 1];
        sigs=[1 1 1 1 1 1];
        % volumetric source value, per zone
        qvf=[1 1 1 0 0 0];
        % incoming flux values
        incf(1:sn) = 0;
        incf((sn/2)+1:sn) = 1;
        % volumetric response value, per zone
        qva=[0 0 1 0 0 0];
        % incoming adj flux values
        inca(1:sn) = 0;
        %Regions to be perturbed. Use value of 1 to specify
        dat.sigaPertRegion=[0 0 0 0 0 0];
        dat.sigsPertRegion=[1 1 1 0 0 0];
        dat.sourcePertRegion=[1 1 1 1 1 1];
        dat.incPertRegion(1:sn) = 1; 
        
case 92 %
        % number of elements per zone
        nel_zone = [ 10 10 10 10 10 10]*4;
        % width of each zone
        width_zone = [1 1 1 1 1 1];
        % sigt/sigs per zone
        sigt=[1 1 1 1 1 1];
        sigs=[1 1 1 1 1 1];
        % volumetric source value, per zone
        qvf=[1 1 1 0 0 0];
        % incoming flux values
        incf(1:sn) = 0;
        incf((sn/2)+1:sn) = 1;
        % volumetric response value, per zone
        qva=[0 0 0 1 0 0];
        % incoming adj flux values
        inca(1:sn) = 0;
        %Regions to be perturbed. Use value of 1 to specify
        dat.sigaPertRegion=[0 0 0 0 0 0];
        dat.sigsPertRegion=[1 1 1 0 0 0];
        dat.sourcePertRegion=[1 1 1 1 1 1];
        dat.incPertRegion(1:sn) = 1;
        
case 93 %Abs to vac
        % number of elements per zone
        nel_zone = [ 10 10 10 10 10 10]*4;
        % width of each zone
        width_zone = [1 1 1 1 1 1];
        % sigt/sigs per zone
        sigt=[1 1 1 1e-8 1e-8 1e-8];
        sigs=[0 0 0 0 0 0];
        % volumetric source value, per zone
        qvf=[1 1 1 0 0 0];
        % incoming flux values
        incf(1:sn) = 0;
        incf((sn/2)+1:sn) = 1;
        % volumetric response value, per zone
        qva=[0 0 1 1 0 0];
        % incoming adj flux values
        inca(1:sn) = 0;
        %Regions to be perturbed. Use value of 1 to specify
        dat.sigaPertRegion=[1 1 1 0 0 0];
        dat.sigsPertRegion=[1 1 1 0 0 0];
        dat.sourcePertRegion=[1 1 1 1 1 1];
        dat.incPertRegion(1:sn) = 1;   
        
case 94 %Abs to vac
        % number of elements per zone
        nel_zone = [ 10 10 10 10 10 10]*4;
        % width of each zone
        width_zone = [1 1 1 1 1 1];
        % sigt/sigs per zone
        sigt=[1 1 1 1e-8 1e-8 1e-8];
        sigs=[0 0 0 0 0 0];
        % volumetric source value, per zone
        qvf=[1 1 1 0 0 0];
        % incoming flux values
        incf(1:sn) = 0;
        incf((sn/2)+1:sn) = 1;
        % volumetric response value, per zone
        qva=[0 0 1 0 0 0];
        % incoming adj flux values
        inca(1:sn) = 0;
        %Regions to be perturbed. Use value of 1 to specify
        dat.sigaPertRegion=[1 1 1 0 0 0];
        dat.sigsPertRegion=[1 1 1 0 0 0];
        dat.sourcePertRegion=[1 1 1 1 1 1];
        dat.incPertRegion(1:sn) = 1;  
      
case 95 %Abs to vac
        % number of elements per zone
        nel_zone = [ 10 10 10 10 10 10]*4;
        % width of each zone
        width_zone = [1 1 1 1 1 1];
        % sigt/sigs per zone
        sigt=[1 1 1 1e-8 1e-8 1e-8];
        sigs=[0 0 0 0 0 0];
        % volumetric source value, per zone
        qvf=[1 1 1 0 0 0];
        % incoming flux values
        incf(1:sn) = 0;
        incf((sn/2)+1:sn) = 1;
        % volumetric response value, per zone
        qva=[0 0 0 1 0 0];
        % incoming adj flux values
        inca(1:sn) = 0;
        %Regions to be perturbed. Use value of 1 to specify
        dat.sigaPertRegion=[1 1 1 0 0 0];
        dat.sigsPertRegion=[1 1 1 0 0 0];
        dat.sourcePertRegion=[1 1 1 1 1 1];
        dat.incPertRegion(1:sn) = 1;  

case 96 %Abs to abs
        % number of elements per zone
        nel_zone = [ 10 10 10 10 10 10]*4;
        % width of each zone
        width_zone = [1 1 1 1 1 1];
        % sigt/sigs per zone
        sigt=[1 1 1 1 1 1];
        sigs=[0 0 0 0 0 0];
        % volumetric source value, per zone
        qvf=[1 1 1 0 0 0];
        % incoming flux values
        incf(1:sn) = 0;
        incf((sn/2)+1:sn) = 1;
        % volumetric response value, per zone
        qva=[0 0 1 1 0 0];
        % incoming adj flux values
        inca(1:sn) = 0;
        %Regions to be perturbed. Use value of 1 to specify
        dat.sigaPertRegion=[1 1 1 0 0 0];
        dat.sigsPertRegion=[1 1 1 0 0 0];
        dat.sourcePertRegion=[1 1 1 1 1 1];
        dat.incPertRegion(1:sn) = 1;
        
        
case 97 %Abs to abs
        % number of elements per zone
        nel_zone = [ 10 10 10 10 10 10]*4;
        % width of each zone
        width_zone = [1 1 1 1 1 1];
        % sigt/sigs per zone
        sigt=[1 1 1 1 1 1];
        sigs=[0 0 0 0 0 0];
        % volumetric source value, per zone
        qvf=[1 1 1 0 0 0];
        % incoming flux values
        incf(1:sn) = 0;
        incf((sn/2)+1:sn) = 1;
        % volumetric response value, per zone
        qva=[0 0 1 0 0 0];
        % incoming adj flux values
        inca(1:sn) = 0;
        %Regions to be perturbed. Use value of 1 to specify
        dat.sigaPertRegion=[1 1 1 0 0 0];
        dat.sigsPertRegion=[1 1 1 0 0 0];
        dat.sourcePertRegion=[1 1 1 1 1 1];
        dat.incPertRegion(1:sn) = 1;
        
case 98 %Abs to abs
        % number of elements per zone
        nel_zone = [ 10 10 10 10 10 10]*4;
        % width of each zone
        width_zone = [1 1 1 1 1 1];
        % sigt/sigs per zone
        sigt=[1 1 1 1 1 1];
        sigs=[0 0 0 0 0 0];
        % volumetric source value, per zone
        qvf=[1 1 1 0 0 0];
        % incoming flux values
        incf(1:sn) = 0;
        incf((sn/2)+1:sn) = 1;
        % volumetric response value, per zone
        qva=[0 0 0 1 0 0];
        % incoming adj flux values
        inca(1:sn) = 0;
        %Regions to be perturbed. Use value of 1 to specify
        dat.sigaPertRegion=[1 1 1 0 0 0];
        dat.sigsPertRegion=[1 1 1 0 0 0];
        dat.sourcePertRegion=[1 1 1 1 1 1];
        dat.incPertRegion(1:sn) = 1;
        
case 99 %Abs to scat
        % number of elements per zone
        nel_zone = [ 10 10 10 10 10 10]*4;
        % width of each zone
        width_zone = [1 1 1 1 1 1];
        % sigt/sigs per zone
        sigt=[1 1 1 1 1 1];
        sigs=[0 0 0 1 1 1];
        % volumetric source value, per zone
        qvf=[1 1 1 0 0 0];
        % incoming flux values
        incf(1:sn) = 0;
        incf((sn/2)+1:sn) = 1;
        % volumetric response value, per zone
        qva=[0 0 1 1 0 0];
        % incoming adj flux values
        inca(1:sn) = 0;
        %Regions to be perturbed. Use value of 1 to specify
        dat.sigaPertRegion=[1 1 1 0 0 0];
        dat.sigsPertRegion=[0 0 0 1 1 1];
        dat.sourcePertRegion=[1 1 1 1 1 1];
        dat.incPertRegion(1:sn) = 1;      
        
case 100 %Abs to scat
        % number of elements per zone
        nel_zone = [ 10 10 10 10 10 10]*4;
        % width of each zone
        width_zone = [1 1 1 1 1 1];
        % sigt/sigs per zone
        sigt=[1 1 1 1 1 1];
        sigs=[0 0 0 1 1 1];
        % volumetric source value, per zone
        qvf=[1 1 1 0 0 0];
        % incoming flux values
        incf(1:sn) = 0;
        incf((sn/2)+1:sn) = 1;
        % volumetric response value, per zone
        qva=[0 0 1 0 0 0];
        % incoming adj flux values
        inca(1:sn) = 0;
        %Regions to be perturbed. Use value of 1 to specify
        dat.sigaPertRegion=[1 1 1 0 0 0];
        dat.sigsPertRegion=[0 0 0 1 1 1];
        dat.sourcePertRegion=[1 1 1 1 1 1];
        dat.incPertRegion(1:sn) = 1;
        
case 101 %Abs to scat
        % number of elements per zone
        nel_zone = [ 10 10 10 10 10 10]*4;
        % width of each zone
        width_zone = [1 1 1 1 1 1];
        % sigt/sigs per zone
        sigt=[1 1 1 1 1 1];
        sigs=[0 0 0 1 1 1];
        % volumetric source value, per zone
        qvf=[1 1 1 0 0 0];
        % incoming flux values
        incf(1:sn) = 0;
        incf((sn/2)+1:sn) = 1;
        % volumetric response value, per zone
        qva=[0 0 0 1 0 0];
        % incoming adj flux values
        inca(1:sn) = 0;
        %Regions to be perturbed. Use value of 1 to specify
        dat.sigaPertRegion=[1 1 1 0 0 0];
        dat.sigsPertRegion=[0 0 0 1 1 1];
        dat.sourcePertRegion=[1 1 1 1 1 1];
        dat.incPertRegion(1:sn) = 1;
        
case 102 %Vac to Abs
        % number of elements per zone
        nel_zone = [ 10 10 10 10 10 10]*4;
        % width of each zone
        width_zone = [1 1 1 1 1 1];
        % sigt/sigs per zone
        sigt=[1e-8 1e-8 1e-8 1 1 1];
        sigs=[0 0 0 0 0 0];
        % volumetric source value, per zone
        qvf=[1 1 1 0 0 0];
        % incoming flux values
        incf(1:sn) = 0;
        incf((sn/2)+1:sn) = 1;
        % volumetric response value, per zone
        qva=[0 0 1 1 0 0];
        % incoming adj flux values
        inca(1:sn) = 0;
        %Regions to be perturbed. Use value of 1 to specify
        dat.sigaPertRegion=[0 0 0 1 1 1];
        dat.sigsPertRegion=[0 0 0 1 1 1];
        dat.sourcePertRegion=[1 1 1 1 1 1];
        dat.incPertRegion(1:sn) = 1;
        
case 103 %Vac to Abs
        % number of elements per zone
        nel_zone = [ 10 10 10 10 10 10]*4;
        % width of each zone
        width_zone = [1 1 1 1 1 1];
        % sigt/sigs per zone
        sigt=[1e-8 1e-8 1e-8 1 1 1];
        sigs=[0 0 0 0 0 0];
        % volumetric source value, per zone
        qvf=[1 1 1 0 0 0];
        % incoming flux values
        incf(1:sn) = 0;
        incf((sn/2)+1:sn) = 1;
        % volumetric response value, per zone
        qva=[0 0 1 0 0 0];
        % incoming adj flux values
        inca(1:sn) = 0;
        %Regions to be perturbed. Use value of 1 to specify
        dat.sigaPertRegion=[0 0 0 1 1 1];
        dat.sigsPertRegion=[0 0 0 1 1 1];
        dat.sourcePertRegion=[1 1 1 1 1 1];
        dat.incPertRegion(1:sn) = 1;
        
case 104 %Vac to Abs
        % number of elements per zone
        nel_zone = [ 10 10 10 10 10 10]*4;
        % width of each zone
        width_zone = [1 1 1 1 1 1];
        % sigt/sigs per zone
        sigt=[1e-8 1e-8 1e-8 1 1 1];
        sigs=[0 0 0 0 0 0];
        % volumetric source value, per zone
        qvf=[1 1 1 0 0 0];
        % incoming flux values
        incf(1:sn) = 0;
        incf((sn/2)+1:sn) = 1;
        % volumetric response value, per zone
        qva=[0 0 0 1 0 0];
        % incoming adj flux values
        inca(1:sn) = 0;
        %Regions to be perturbed. Use value of 1 to specify
        dat.sigaPertRegion=[0 0 0 1 1 1];
        dat.sigsPertRegion=[0 0 0 1 1 1];
        dat.sourcePertRegion=[1 1 1 1 1 1];
        dat.incPertRegion(1:sn) = 1;
        
case 105 %Vac to Scat
        % number of elements per zone
        nel_zone = [ 10 10 10 10 10 10]*4;
        % width of each zone
        width_zone = [1 1 1 1 1 1];
        % sigt/sigs per zone
        sigt=[1e-8 1e-8 1e-8 1 1 1];
        sigs=[0 0 0 1 1 1];
        % volumetric source value, per zone
        qvf=[1 1 1 0 0 0];
        % incoming flux values
        incf(1:sn) = 0;
        incf((sn/2)+1:sn) = 1;
        % volumetric response value, per zone
        qva=[0 0 1 1 0 0];
        % incoming adj flux values
        inca(1:sn) = 0;
        %Regions to be perturbed. Use value of 1 to specify
        dat.sigaPertRegion=[0 0 0 1 1 1];
        dat.sigsPertRegion=[0 0 0 1 1 1];
        dat.sourcePertRegion=[1 1 1 1 1 1];
        dat.incPertRegion(1:sn) = 1;
        
case 106 %Vac to Scat
        % number of elements per zone
        nel_zone = [ 10 10 10 10 10 10]*4;
        % width of each zone
        width_zone = [1 1 1 1 1 1];
        % sigt/sigs per zone
        sigt=[1e-8 1e-8 1e-8 1 1 1];
        sigs=[0 0 0 1 1 1];
        % volumetric source value, per zone
        qvf=[1 1 1 0 0 0];
        % incoming flux values
        incf(1:sn) = 0;
        incf((sn/2)+1:sn) = 1;
        % volumetric response value, per zone
        qva=[0 0 1 0 0 0];
        % incoming adj flux values
        inca(1:sn) = 0;
        %Regions to be perturbed. Use value of 1 to specify
        dat.sigaPertRegion=[0 0 0 1 1 1];
        dat.sigsPertRegion=[0 0 0 1 1 1];
        dat.sourcePertRegion=[1 1 1 1 1 1];
        dat.incPertRegion(1:sn) = 1;
        
case 107 %Vac to Scat
        % number of elements per zone
        nel_zone = [ 10 10 10 10 10 10]*4;
        % width of each zone
        width_zone = [1 1 1 1 1 1];
        % sigt/sigs per zone
        sigt=[1e-8 1e-8 1e-8 1 1 1];
        sigs=[0 0 0 1 1 1];
        % volumetric source value, per zone
        qvf=[1 1 1 0 0 0];
        % incoming flux values
        incf(1:sn) = 0;
        incf((sn/2)+1:sn) = 1;
        % volumetric response value, per zone
        qva=[0 0 0 1 0 0];
        % incoming adj flux values
        inca(1:sn) = 0;
        %Regions to be perturbed. Use value of 1 to specify
        dat.sigaPertRegion=[0 0 0 1 1 1];
        dat.sigsPertRegion=[0 0 0 1 1 1];
        dat.sourcePertRegion=[1 1 1 1 1 1];
        dat.incPertRegion(1:sn) = 1;        
        
%%%%Testing something with source additivity
case 108 %
        % number of elements per zone
        nel_zone = [ 10 10 10 10 10 10]*4;
        % width of each zone
        width_zone = [1 1 1 1 1 1];
        % sigt/sigs per zone
        sigt=[1 1 1 1 1 1];
        sigs=[1 1 1 0 0 0];
        % volumetric source value, per zone
        qvf=[1 1 1 0 0 0];
        % incoming flux values
        incf(1:sn) = 0;
        incf((sn/2)+1:sn) = 0;
        % volumetric response value, per zone
        qva=[0 0 1 1 0 0];
        % incoming adj flux values
        inca(1:sn) = 0;
        %Regions to be perturbed. Use value of 1 to specify
        dat.sigaPertRegion=[0 0 0 1 1 1];
        dat.sigsPertRegion=[1 1 1 0 0 0];
        dat.sourcePertRegion=[1 1 1 1 1 1];
        dat.incPertRegion(1:sn) = 1; 

case 109 %
        % number of elements per zone
        nel_zone = [ 10 10 10 10 10 10]*4;
        % width of each zone
        width_zone = [1 1 1 1 1 1];
        % sigt/sigs per zone
        sigt=[1 1 1 1 1 1];
        sigs=[1 1 1 0 0 0];
        % volumetric source value, per zone
        qvf=[1 1 1 0 0 0];
        % incoming flux values
        incf(1:sn) = 0;
        incf((sn/2)+1:sn) = 1;
        % volumetric response value, per zone
        qva=[0 0 1 1 0 0];
        % incoming adj flux values
        inca(1:sn) = 0;
        %Regions to be perturbed. Use value of 1 to specify
        dat.sigaPertRegion=[0 0 0 1 1 1];
        dat.sigsPertRegion=[1 1 1 0 0 0];
        dat.sourcePertRegion=[1 1 1 1 1 1];
        dat.incPertRegion(1:sn) = 1; 
        
   case 110 %Uniform, 0 inc flux. Response in middle
        % number of elements per zone
        nel_zone = [ 10 10 10 10 10]*4;
        % width of each zone
        width_zone = [ 2 2 2 2 2];
        % sigt/sigs per zone
        sigt=[2 2 2 2 2];
        sigs=[1 1 1 1 1];
        % volumetric source value, per zone
        qvf=[1 1 1 1 1];
        % incoming flux values
        incf(1:sn) = 0;
        % volumetric response value, per zone
        qva=[0 0 1 0 0];
        % incoming adj flux values
        inca(1:sn) = 0;
        %Regions to be perturbed. Use value of 1 to specify
        dat.sigaPertRegion=[1 0 1 0 1];
        dat.sigsPertRegion=[1 0 1 0 1];
        dat.sourcePertRegion=[1 0 1 0 1];
        dat.incPertRegion(1:sn)=1;
        
   case 111 %Uniform, 0 inc flux. Response in middle
        % number of elements per zone
        nel_zone = [ 10 10 10 10 10]*4;
        % width of each zone
        width_zone = [ 2 2 2 2 2];
        % sigt/sigs per zone
        sigt=[2 2 2 2 2];
        sigs=[1 1 1 1 1];
        % volumetric source value, per zone
        qvf=[1 1 1 1 1];
        % incoming flux values
        incf(1:sn) = 0;
        % volumetric response value, per zone
        qva=[0 0 1 0 0];
        % incoming adj flux values
        inca(1:sn) = 0;
        %Regions to be perturbed. Use value of 1 to specify
        dat.sigaPertRegion=[0 0 1 1 1];
        dat.sigsPertRegion=[0 0 1 1 1];
        dat.sourcePertRegion=[0 0 1 1 1];
        dat.incPertRegion(1:sn)=1;

   case 112 %Uniform, 0 inc flux. Response in middle
        % number of elements per zone
        nel_zone = [ 10 10 10 10 10]*4;
        % width of each zone
        width_zone = [ 2 2 2 2 2];
        % sigt/sigs per zone
        sigt=[2 2 2 2 2];
        sigs=[1 1 1 1 1];
        % volumetric source value, per zone
        qvf=[1 1 1 1 1];
        % incoming flux values
        incf(1:sn) = 0;
        % volumetric response value, per zone
        qva=[0 0 1 0 0];
        % incoming adj flux values
        inca(1:sn) = 0;
        %Regions to be perturbed. Use value of 1 to specify
        dat.sigaPertRegion=[1 0 -1 0 1];
        dat.sigsPertRegion=[1 0 -1 0 1];
        dat.sourcePertRegion=[1 0 -1 0 1];
        dat.incPertRegion(1:sn)=1;
        
    otherwise
        error('problem ID %g is unknown',pb_ID);
end
return
end