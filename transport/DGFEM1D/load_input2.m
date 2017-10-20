function [nel_zone,width_zone,sigt,sigs,qvf,incf,qva,inca]=load_input2(pb_ID)
% load XS (total and isotropic scattering) for each material
% and load each volumetric external (isotropic) source intensity
% non-isotropic data is not considered here

global dat snq npar

sn=snq.n_dir;
mu=snq.mu;

%%%%%%%%%%%%%%%%%%%%%% select input data

switch pb_ID
case 60 %
        % number of elements per zone
        nel_zone = [ 10 10 10 10 10 10]*4;
        % width of each zone
        width_zone = [1 1 1 1 1 1];
        % sigt/sigs per zone
        sigt=[1 1 1 1e-8 1e-8 1e-8];
        sigs=[1 1 1 0 0 0];
        % volumetric source value, per zone
        qvf=[0 0 0 0 0 0];
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
        dat.sourcePertRegion=[0 0 0 0 0 0];
        dat.incPertRegion(1:sn) = 1;
        
case 61 %
        % number of elements per zone
        nel_zone = [ 10 10 10 10 10 10]*4;
        % width of each zone
        width_zone = [1 1 1 1 1 1];
        % sigt/sigs per zone
        sigt=[1 1 1 1e-8 1e-8 1e-8];
        sigs=[1 1 1 0 0 0];
        % volumetric source value, per zone
        qvf=[0 0 0 0 0 0];
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
        dat.sourcePertRegion=[0 0 0 0 0 0];
        dat.incPertRegion(1:sn) = 1;
       
case 62 %
        % number of elements per zone
        nel_zone = [ 10 10 10 10 10 10]*4;
        % width of each zone
        width_zone = [1 1 1 1 1 1];
        % sigt/sigs per zone
        sigt=[1 1 1 1e-8 1e-8 1e-8];
        sigs=[1 1 1 0 0 0];
        % volumetric source value, per zone
        qvf=[0 0 0 0 0 0];
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
        dat.sourcePertRegion=[0 0 0 0 0 0];
        dat.incPertRegion(1:sn) = 1;   
        
case 63 %
        % number of elements per zone
        nel_zone = [ 10 10 10 10 10 10]*4;
        % width of each zone
        width_zone = [1 1 1 1 1 1];
        % sigt/sigs per zone
        sigt=[1 1 1 1 1 1];
        sigs=[1 1 1 0 0 0];
        % volumetric source value, per zone
        qvf=[0 0 0 0 0 0];
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
        dat.sourcePertRegion=[0 0 0 0 0 0];
        dat.incPertRegion(1:sn) = 1; 
        
case 64 %
        % number of elements per zone
        nel_zone = [ 10 10 10 10 10 10]*4;
        % width of each zone
        width_zone = [1 1 1 1 1 1];
        % sigt/sigs per zone
        sigt=[1 1 1 1 1 1];
        sigs=[1 1 1 0 0 0];
        % volumetric source value, per zone
        qvf=[0 0 0 0 0 0];
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
        dat.sourcePertRegion=[0 0 0 0 0 0];
        dat.incPertRegion(1:sn) = 1; 
        
case 65 %
        % number of elements per zone
        nel_zone = [ 10 10 10 10 10 10]*4;
        % width of each zone
        width_zone = [1 1 1 1 1 1];
        % sigt/sigs per zone
        sigt=[1 1 1 1 1 1];
        sigs=[1 1 1 0 0 0];
        % volumetric source value, per zone
        qvf=[0 0 0 0 0 0];
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
        dat.sourcePertRegion=[0 0 0 0 0 0];
        dat.incPertRegion(1:sn) = 1;
        
case 66 %
        % number of elements per zone
        nel_zone = [ 10 10 10 10 10 10]*4;
        % width of each zone
        width_zone = [1 1 1 1 1 1];
        % sigt/sigs per zone
        sigt=[1 1 1 1 1 1];
        sigs=[1 1 1 1 1 1];
        % volumetric source value, per zone
        qvf=[0 0 0 0 0 0];
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
        dat.sourcePertRegion=[0 0 0 0 0 0];
        dat.incPertRegion(1:sn) = 1;
        
case 67 %
        % number of elements per zone
        nel_zone = [ 10 10 10 10 10 10]*4;
        % width of each zone
        width_zone = [1 1 1 1 1 1];
        % sigt/sigs per zone
        sigt=[1 1 1 1 1 1];
        sigs=[1 1 1 1 1 1];
        % volumetric source value, per zone
        qvf=[0 0 0 0 0 0];
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
        dat.sourcePertRegion=[0 0 0 0 0 0];
        dat.incPertRegion(1:sn) = 1; 
        
case 68 %
        % number of elements per zone
        nel_zone = [ 10 10 10 10 10 10]*4;
        % width of each zone
        width_zone = [1 1 1 1 1 1];
        % sigt/sigs per zone
        sigt=[1 1 1 1 1 1];
        sigs=[1 1 1 1 1 1];
        % volumetric source value, per zone
        qvf=[0 0 0 0 0 0];
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
        dat.sourcePertRegion=[0 0 0 0 0 0];
        dat.incPertRegion(1:sn) = 1;
        
case 69 %Abs to vac
        % number of elements per zone
        nel_zone = [ 10 10 10 10 10 10]*4;
        % width of each zone
        width_zone = [1 1 1 1 1 1];
        % sigt/sigs per zone
        sigt=[1 1 1 1e-8 1e-8 1e-8];
        sigs=[0 0 0 0 0 0];
        % volumetric source value, per zone
        qvf=[0 0 0 0 0 0];
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
        dat.sourcePertRegion=[0 0 0 0 0 0];
        dat.incPertRegion(1:sn) = 1;   
        
case 70 %Abs to vac
        % number of elements per zone
        nel_zone = [ 10 10 10 10 10 10]*4;
        % width of each zone
        width_zone = [1 1 1 1 1 1];
        % sigt/sigs per zone
        sigt=[1 1 1 1e-8 1e-8 1e-8];
        sigs=[0 0 0 0 0 0];
        % volumetric source value, per zone
        qvf=[0 0 0 0 0 0];
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
        dat.sourcePertRegion=[0 0 0 0 0 0];
        dat.incPertRegion(1:sn) = 1;  
      
case 71 %Abs to vac
        % number of elements per zone
        nel_zone = [ 10 10 10 10 10 10]*4;
        % width of each zone
        width_zone = [1 1 1 1 1 1];
        % sigt/sigs per zone
        sigt=[1 1 1 1e-8 1e-8 1e-8];
        sigs=[0 0 0 0 0 0];
        % volumetric source value, per zone
        qvf=[0 0 0 0 0 0];
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
        dat.sourcePertRegion=[0 0 0 0 0 0];
        dat.incPertRegion(1:sn) = 1;  

case 72 %Abs to abs
        % number of elements per zone
        nel_zone = [ 10 10 10 10 10 10]*4;
        % width of each zone
        width_zone = [1 1 1 1 1 1];
        % sigt/sigs per zone
        sigt=[1 1 1 1 1 1];
        sigs=[0 0 0 0 0 0];
        % volumetric source value, per zone
        qvf=[0 0 0 0 0 0];
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
        dat.sourcePertRegion=[0 0 0 0 0 0];
        dat.incPertRegion(1:sn) = 1;
        
        
case 73 %Abs to abs
        % number of elements per zone
        nel_zone = [ 10 10 10 10 10 10]*4;
        % width of each zone
        width_zone = [1 1 1 1 1 1];
        % sigt/sigs per zone
        sigt=[1 1 1 1 1 1];
        sigs=[0 0 0 0 0 0];
        % volumetric source value, per zone
        qvf=[0 0 0 0 0 0];
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
        dat.sourcePertRegion=[0 0 0 0 0 0];
        dat.incPertRegion(1:sn) = 1;
        
case 74 %Abs to abs
        % number of elements per zone
        nel_zone = [ 10 10 10 10 10 10]*4;
        % width of each zone
        width_zone = [1 1 1 1 1 1];
        % sigt/sigs per zone
        sigt=[1 1 1 1 1 1];
        sigs=[0 0 0 0 0 0];
        % volumetric source value, per zone
        qvf=[0 0 0 0 0 0];
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
        dat.sourcePertRegion=[0 0 0 0 0 0];
        dat.incPertRegion(1:sn) = 1;
        
case 75 %Abs to scat
        % number of elements per zone
        nel_zone = [ 10 10 10 10 10 10]*4;
        % width of each zone
        width_zone = [1 1 1 1 1 1];
        % sigt/sigs per zone
        sigt=[1 1 1 1 1 1];
        sigs=[0 0 0 1 1 1];
        % volumetric source value, per zone
        qvf=[0 0 0 0 0 0];
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
        dat.sourcePertRegion=[0 0 0 0 0 0];
        dat.incPertRegion(1:sn) = 1;      
        
case 76 %Abs to scat
        % number of elements per zone
        nel_zone = [ 10 10 10 10 10 10]*4;
        % width of each zone
        width_zone = [1 1 1 1 1 1];
        % sigt/sigs per zone
        sigt=[1 1 1 1 1 1];
        sigs=[0 0 0 1 1 1];
        % volumetric source value, per zone
        qvf=[0 0 0 0 0 0];
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
        dat.sourcePertRegion=[0 0 0 0 0 0];
        dat.incPertRegion(1:sn) = 1;
        
case 77 %Abs to scat
        % number of elements per zone
        nel_zone = [ 10 10 10 10 10 10]*4;
        % width of each zone
        width_zone = [1 1 1 1 1 1];
        % sigt/sigs per zone
        sigt=[1 1 1 1 1 1];
        sigs=[0 0 0 1 1 1];
        % volumetric source value, per zone
        qvf=[0 0 0 0 0 0];
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
        dat.sourcePertRegion=[0 0 0 0 0 0];
        dat.incPertRegion(1:sn) = 1;
        
case 78 %Vac to Abs
        % number of elements per zone
        nel_zone = [ 10 10 10 10 10 10]*4;
        % width of each zone
        width_zone = [1 1 1 1 1 1];
        % sigt/sigs per zone
        sigt=[1e-8 1e-8 1e-8 1 1 1];
        sigs=[0 0 0 0 0 0];
        % volumetric source value, per zone
        qvf=[0 0 0 0 0 0];
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
        dat.sourcePertRegion=[0 0 0 0 0 0];
        dat.incPertRegion(1:sn) = 1;
        
case 79 %Vac to Abs
        % number of elements per zone
        nel_zone = [ 10 10 10 10 10 10]*4;
        % width of each zone
        width_zone = [1 1 1 1 1 1];
        % sigt/sigs per zone
        sigt=[1e-8 1e-8 1e-8 1 1 1];
        sigs=[0 0 0 0 0 0];
        % volumetric source value, per zone
        qvf=[0 0 0 0 0 0];
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
        dat.sourcePertRegion=[0 0 0 0 0 0];
        dat.incPertRegion(1:sn) = 1;
        
case 80 %Vac to Abs
        % number of elements per zone
        nel_zone = [ 10 10 10 10 10 10]*4;
        % width of each zone
        width_zone = [1 1 1 1 1 1];
        % sigt/sigs per zone
        sigt=[1e-8 1e-8 1e-8 1 1 1];
        sigs=[0 0 0 0 0 0];
        % volumetric source value, per zone
        qvf=[0 0 0 0 0 0];
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
        dat.sourcePertRegion=[0 0 0 0 0 0];
        dat.incPertRegion(1:sn) = 1;
        
case 81 %Vac to Scat
        % number of elements per zone
        nel_zone = [ 10 10 10 10 10 10]*4;
        % width of each zone
        width_zone = [1 1 1 1 1 1];
        % sigt/sigs per zone
        sigt=[1e-8 1e-8 1e-8 1 1 1];
        sigs=[0 0 0 1 1 1];
        % volumetric source value, per zone
        qvf=[0 0 0 0 0 0];
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
        dat.sourcePertRegion=[0 0 0 0 0 0];
        dat.incPertRegion(1:sn) = 1;
        
case 82 %Vac to Scat
        % number of elements per zone
        nel_zone = [ 10 10 10 10 10 10]*4;
        % width of each zone
        width_zone = [1 1 1 1 1 1];
        % sigt/sigs per zone
        sigt=[1e-8 1e-8 1e-8 1 1 1];
        sigs=[0 0 0 1 1 1];
        % volumetric source value, per zone
        qvf=[0 0 0 0 0 0];
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
        dat.sourcePertRegion=[0 0 0 0 0 0];
        dat.incPertRegion(1:sn) = 1;
        
case 83 %Vac to Scat
        % number of elements per zone
        nel_zone = [ 10 10 10 10 10 10]*4;
        % width of each zone
        width_zone = [1 1 1 1 1 1];
        % sigt/sigs per zone
        sigt=[1e-8 1e-8 1e-8 1 1 1];
        sigs=[0 0 0 1 1 1];
        % volumetric source value, per zone
        qvf=[0 0 0 0 0 0];
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
        dat.sourcePertRegion=[0 0 0 0 0 0];
        dat.incPertRegion(1:sn) = 1;        

        
    otherwise
        error('problem ID %g is unknown',pb_ID);
end
return
end