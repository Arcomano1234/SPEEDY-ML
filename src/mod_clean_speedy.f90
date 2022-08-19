module mod_clean_speedy
   !We need to flush the memeory of speedy 
   !This module sets every speedy variable to zero

   use mod_cli_land
   use mod_cli_sea 
   use mod_cplcon_sea
   use mod_cpl_land_model
   use mod_cplvar_sea
   use mod_dyncon1
 
   subroutine clean_up_speedy() 
      !mod_cli_land vars
      fmask_l = 0.0
      bmask_l = 0.0
      stl12 = 0.0
      snowd12 = 0.0
      soilw12 = 0.0

      !mod_cli_sea.f90 vars
      fmask_s = 0.0
      bmask_s = 0.0
      deglat_s = 0.0
      sst12 = 0.0
 
      !mod_cplcon_sea vars
      rhcaps = 0.0
      rhcapi = 0.0 
      cdsea = 0.0 
      cdice = 0.0 
    
      !mod_cpl_land_model
      rhcapl = 0.0
      cdland = 0.0 
      vland_input = 0.0 
      vland_output = 0.0

      !mod_cplvar_sea
      vsea_input = 0.0 
      vsea_output = 0.0

      !mod_dyncon1
      hsg = 0.0
      dhs = 0.0 
      fsg = 0.0
      dhsr = 0.0
      fsgr = 0.0
      radang = 0.0
      gsin = 0.0 
      gcos = 0.0
      coriol = 0.0
      xgeop1 = 0.0
      xgeop2 = 0.0
   

      

    
