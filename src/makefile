COMPILE1=mpif90
COMPILE2=mpiifort

COMLIB1=

COMOTT2=  -traceback -mkl -lnetcdff -larpack -heap-arrays 10000 -mcmodel=medium -nowarn -O3 -lhdf5_fortran -lhdf5 -r8 -msse3 -traceback -standard-realloc-lhs -ipo -fpp -lmkl_intel_lp64 -lmkl_intel_thread -liomp5 -lpthread #-g -check all -check noarg_temp_created -warn -traceback -debug extended -mkl -larpack -heap-arrays 10000 -mcmodel=medium -lhdf5_fortran -lhdf5 -r8 -msse3 -traceback -standard-realloc-lhs -ipo -fpp -lmkl_intel_lp64 -lmkl_intel_thread -liomp5 -lpthread -I/sw/eb/sw/netCDF-Fortran/4.5.2-iimpi-2020a/include -L/sw/eb/sw/netCDF-Fortran/4.5.2-iimpi-2020a/lib -lnetcdff -L/sw/eb/sw/netCDF/4.7.4-iimpi-2020a/lib64 -L/sw/eb/sw/iccifort/2020.1.217/lib/intel64 -lnetcdf -lnetcdf -lm -liomp5 -lpthread #-g -check all -check noarg_temp_created -warn -traceback -debug extended -mkl -lnetcdff -larpack -heap-arrays 10000 -mcmodel=medium -O3 -lhdf5_fortran -lhdf5 -r8 -msse3 -traceback -standard-realloc-lhs -ipo -fpp -lmkl_intel_lp64 -lmkl_intel_thread -liomp5 -lpthread #-axCORE-AVX2\\

DEBUG_FLAGS= -g -fpe0 -check bounds -warn -traceback -debug extended -mkl -lnetcdff -larpack -heap-arrays 10000 -mcmodel=medium -lhdf5_fortran -lhdf5 -r8 -lmkl_intel_lp64 -lmkl_intel_thread -liomp5 -lpthread
#COMOTT2= -r8 -traceback 
#-O3 -ffast-math -finline-matmul-limit=0 
#-g -Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace -fcheck=all -ffpe-trap=zero,overflow,underflow -fmax-stack-var-size=65535
COMOTT1 = -O3 -ffast-math -finline-matmul-limit=0 -finit-real=zero -finit-logical=false -finit-integer=zero -mcmodel=medium -ffree-line-length-none -fallow-argument-mismatch -fconvert=swap -fdefault-real-8 -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl  -m64  -I"${MKLROOT}/include" -I${EBROOTNETCDFMINFORTRAN}/include -L${EBROOTNETCDFMINFORTRAN}/lib64 -lnetcdff -L${EBROOTNETCDF}/lib64 -L${EBROOTGCCCORE}/lib64 -L${EBROOTGCCCORE}/lib -lnetcdf -lnetcdf -ldl -lm -lm -lpthread -I${HDF5_DIR}/include -L${HDF5_DIR}/lib -lhdf5_fortran -lhdf5 -larpack #g  -Wall -fcheck=al -O3 -ffast-math -march=native -funroll-loops -fcheck=all -g -fbacktrace -O3 -ffast-math -funroll-loops
#COMOTT1= -O3 -ffast-math -fdefault-integer-8 -ffree-line-length-512 -fallow-argument-mismatch -fbacktrace -fconvert=swap -fdefault-real-8 -larpack -DMKL_ILP64 -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_gf_ilp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl -m64  -I"${MKLROOT}/include" -I${EBROOTNETCDFMINFORTRAN}/include -L${EBROOTNETCDFMINFORTRAN}/lib -lnetcdff -L${EBROOTNETCDF}/lib64 -L${EBROOTGCCCORE}/lib64 -L${EBROOTGCCCORE}/lib -lnetcdf -lnetcdf -ldl -lm -lm -lpthread -I${HDF5_DIR}/include -L${HDF5_DIR}/lib -lhdf5_fortran -lhdf5 std=f2008 -O3 -ffast-math -march=native -funroll-loops

#COMOTT1= -O3 -ffast-math -fdefault-integer-8 -fallow-argument-mismatch -fbacktrace -fconvert=swap -fdefault-real-8 -DMKL_ILP64 -L/sw/eb/sw/arpack-ng/3.8.0-foss-2021a/lib64 -larpack -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_gf_ilp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl  -m64  -I"${MKLROOT}/include" -liomp5 -lpthread -lm -lmkl_gf_ilp64 -lmkl_intel_thread -lmkl_core -ffree-line-length-512 -I${EBROOTNETCDFMINFORTRAN}/include -L${EBROOTNETCDFMINFORTRAN}/lib -lnetcdff -L${EBROOTNETCDF}/lib64 -L${EBROOTGCCCORE}/lib64 -L${EBROOTGCCCORE}/lib -lnetcdf -lnetcdf -ldl -lm -lm -lpthread -I${HDF5_DIR}/include -L${HDF5_DIR}/lib -lhdf5_fortran -lhdf5#-fcheck=all -fbacktrace -fconvert=swap -fdefault-real-8 -larpack -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_gf_ilp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl  -m64  -I"${MKLROOT}/include" -liomp5 -lpthread -lm -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -ffree-line-length-512 -I${EBROOTNETCDFMINFORTRAN}/include -L${EBROOTNETCDFMINFORTRAN}/lib -lnetcdff -L${EBROOTNETCDF}/lib64 -L${EBROOTGCCCORE}/lib64 -L${EBROOTGCCCORE}/lib -lnetcdf -lnetcdf -ldl -lm -lm -lpthread -I${HDF5_DIR}/include -L${HDF5_DIR}/lib -lhdf5_fortran -lhdf5 -O3 -ffast-math -fallow-argument-mismatch -Wall -fcheck=all

 COMCONV2= #-convert big_endian
 COMCONV1= -fconvert=swap


ifeq ($(COMPILER),gcc2021)
        COMPILE=$(COMPILE1)
        COMOTT=$(COMOTT1) 
        COMCONV=$(COMCONV1) 
        COMLIB=$(COMLIB1)

        MKL_SPARSE_FILE=mkl_spblas.mod_gcc2021
else ifeq ($(COMPILER),gcc2020)
        COMPILE=$(COMPILE1)
        COMOTT=$(COMOTT1)
        COMCONV=$(COMCONV1)
        COMLIB=$(COMLIB1)

        MKL_SPARSE_FILE=mkl_spblas.mod_gcc2020

else ifeq ($(COMPILER),gcc2018)
        COMPILE=$(COMPILE1)
        COMOTT=$(COMOTT1)
        COMCONV=$(COMCONV1)
        COMLIB=$(COMLIB1)

        MKL_SPARSE_FILE=mkl_spblas.mod_gcc2018

else ifeq ($(COMPILER),intel2021)
        COMPILE=$(COMPILE2)
        COMOTT=$(COMOTT2)
        COMCONV=$(COMCONV2)
        COMLIB=$(COMLIB2)

        MKL_SPARSE_FILE=mkl_spblas.mod_intel_2021

else
        COMPILE=$(COMPILE2)
        COMOTT=$(COMOTT2)
        COMCONV=$(COMCONV2)
        COMLIB=$(COMLIB2)

        MKL_SPARSE_FILE=mkl_spblas.mod_intel_2020
endif

FILES= \
       at_gcm.o \
       cpl_bcinterp.o \
       cpl_land.o \
       cpl_main_interface.o \
       cpl_sea.o \
       cpl_sea_model.o \
       dyn_geop.o \
       dyn_grtend.o \
       dyn_implic.o \
       dyn_sptend.o \
       dyn_step.o \
       dyn_stloop.o \
       ini_agcm_init.o \
       ini_impint.o \
       ini_indyns.o \
       ini_inbcon.o \
       ini_iniatm.o \
       ini_inphys.o \
       ini_inirdf.o \
       ini_invars.o \
       ini_stepone.o \
       ini_fordate.o \
	   mod_atparam.o \
	   mod_cli_land.o \
	   mod_cli_sea.o \
	   mod_cnvcon.o \
	   mod_cplcon_sea.o \
	   mod_cplvar_sea.o \
	   mod_cpl_flags.o \
       mod_cpl_land_model.o \
	   mod_date.o \
	   mod_dyncon0.o \
	   mod_dyncon1.o \
	   mod_dyncon2.o \
	   mod_dynvar.o \
	   mod_fft.o \
	   mod_flx_land.o \
	   mod_flx_sea.o \
	   mod_hdifcon.o \
	   mod_lflags.o \
	   mod_lsccon.o \
	   mod_physcon.o \
	   mod_physvar.o \
	   mod_radcon.o \
	   mod_randfor.o \
	   mod_sflcon.o \
	   mod_spectral.o \
	   mod_sppt.o \
	   mod_surfcon.o \
	   mod_tmean.o \
	   mod_tsteps.o \
	   mod_var_land.o \
	   mod_var_sea.o \
	   mod_vdicon.o \
       phy_convmf.o \
       phy_lscond.o \
       phy_phypar.o \
       phy_radiat.o \
       phy_shtorh.o \
       phy_suflux.o \
       phy_vdifsc.o \
       ppo_diagns.o \
       ppo_iogrid.o \
       ppo_restart.o \
       ppo_setctl.o \
       ppo_setgrd.o \
       ppo_tminc.o \
       ppo_dmflux.o \
       ppo_tmout.o \
       ppo_dmout.o \
       spe_matinv.o \
       spe_spectral.o \
       spe_subfft_fftpack.o \
       stringtype.o \
       mod_utilities.o \
       mod_calendar.o \
       mod_linalg.o \
       res_domain.o \
       mpires.o \
       mod_reservoir.o \
       mod_slab_ocean_reservoir.o \
       mod_io.o \
       speedy_res_interface.o \
       parallelmain.o

.SUFFIXES: .f .f90 .o

%.o: %.mod

.f90.o:
	$(COMPILE) $(COMCONV) $(COMOTT) -c $<

.f.o:
	$(COMPILE) $(COMCONV) $(COMOTT) -c $<

imp.exe: move_mod $(FILES)  
	$(COMPILE) $(COMCONV) $(COMOTT) $(FILES) -o imp.exe $(COMLIB)

imp_gcc.exe: $(FILES)
	$(COMPILE) $(COMCONV) $(COMOTT) $(FILES) -o imp_gcc.exe $(COMLIB)

test.exe: $(FILES)
	$(COMPILE) $(COMCONV) $(COMOTT) $(FILES) -o test.exe $(COMLIB)

hybrid.exe: $(FILES)
	$(COMPILE) $(COMCONV) $(COMOTT) $(FILES) -o hybrid.exe $(COMLIB)

debug: $(FILES)	
	$(COMPILE) $(COMCONV) $(COMOTT) $(FILES) -o imp.exe $(COMLIB)

.PHONY: clean
clean:
	rm -f *.o 
	find *.mod -not -name mkl_spblas.mod -delete

move_mod: 
	cp $(MKL_SPARSE_FILE) mkl_spblas.mod

mod_utilities.o        : mod_utilities.f90
stringtype.o           : stringtype.f90
mod_calendar.o         : mod_calendar.f90 mod_utilities.o
mod_linalg.o           : mod_linalg.f90 mod_utilities.o
res_domain.o           : res_domain.f90 mod_utilities.o
mod_io.o               : mod_io.f90 mod_utilities.o res_domain.o stringtype.o
mpires.o               : mpires.f90 mod_utilities.o res_domain.o mod_io.o stringtype.o mod_calendar.o
mod_reservoir.o        : mod_reservoir.f90 mod_utilities.o res_domain.o mod_linalg.o mpires.o speedy_res_interface.o
mod_slab_ocean_reservoir.o : mod_utilities.o res_domain.o mod_linalg.o mpires.o speedy_res_interface.o
speedy_res_interface.o : speedy_res_interface.f90 mod_utilities.o mod_atparam.o mod_physvar.o \
                           mod_tsteps.o mod_calendar.o  mpires.o mod_io.o res_domain.o \
                           at_gcm.o
parallelmain.o         : parallelmain.f90 mpires.o mod_reservoir.o mod_utilities.o \
                            speedy_res_interface.o res_domain.o 


at_gcm.o               : at_gcm.f90 mod_tsteps.o mod_date.o mod_utilities.o
#addflx.o               : addflx.f
#atm2land.o             : atm2land.f
#atm2sea.o              : atm2sea.f
#get_from_coupler.o     : get_from_coupler.f
cpl_bcinterp.o         : cpl_bcinterp.f90
cpl_land.o             : cpl_land.f90 mod_cpl_flags.o mod_atparam.o mod_cpl_land_model.o mod_flx_land.o\
						 mod_cli_land.o mod_date.o mod_var_land.o
mod_cpl_land_model.o   : mod_cpl_land_model.f90 mod_atparam.o
cpl_main_interface.o   : cpl_main_interface.f90 mod_atparam.o mod_cpl_land_model.o mod_surfcon.o\
	                     mod_cli_land.o mod_cli_sea.o
cpl_sea.o              : cpl_sea.f90 mod_cpl_flags.o mod_atparam.o mod_cplvar_sea.o mod_date.o mod_flx_sea.o\
	                     mod_cli_sea.o mod_var_sea.o mpires.o
cpl_sea_model.o        : cpl_sea_model.f90 mod_atparam.o mod_cplcon_sea.o mod_cplvar_sea.o
dyn_geop.o             : dyn_geop.f90 mod_atparam.o mod_dynvar.o mod_dyncon1.o
dyn_grtend.o           : dyn_grtend.f90 mod_atparam.o mod_dynvar.o mod_dyncon1.o mod_dyncon2.o
dyn_implic.o           : dyn_implic.f90 mod_atparam.o mod_dyncon1.o mod_dyncon1.o mod_dyncon2.o
dyn_sptend.o           : dyn_sptend.f90 mod_atparam.o mod_dynvar.o mod_dyncon1.o mod_dyncon2.o 
dyn_step.o             : dyn_step.f90 mod_tsteps.o mod_dyncon0.o mod_atparam.o mod_dynvar.o mod_hdifcon.o mod_reservoir.o
dyn_stloop.o           : dyn_stloop.f90 mod_lflags.o mod_tsteps.o mod_date.o speedy_res_interface.o
#ini_coupler.o          : ini_coupler.f
ini_agcm_init.o        : ini_agcm_init.f90 mod_cpl_flags.o mod_tsteps.o mod_atparam.o mod_date.o
ini_impint.o           : ini_impint.f90 mod_dyncon0.o mod_atparam.o mod_dyncon1.o mod_dyncon2.o mod_hdifcon.o
ini_indyns.o           : ini_indyns.f90 mod_tsteps.o mod_dyncon0.o mod_atparam.o mod_dyncon1.o mod_hdifcon.o\
	                     mod_spectral.o
ini_inbcon.o           : ini_inbcon.f90 mod_cpl_flags.o mod_tsteps.o mod_dyncon0.o mod_atparam.o mod_dyncon1.o\
						 mod_surfcon.o mod_cli_land.o mod_cli_sea.o 
ini_iniatm.o           : ini_iniatm.f90 mod_tsteps.o mod_atparam.o mod_dyncon1.o mod_tmean.o mod_date.o
ini_inphys.o           : ini_inphys.f90 mod_atparam.o mod_physcon.o
ini_inirdf.o           : ini_inirdf.f90 mod_atparam.o mod_physcon.o mod_randfor.o
ini_invars.o           : ini_invars.f90 mod_tsteps.o mod_dyncon0.o mod_atparam.o mod_dynvar.o mod_dyncon1.o\
						 mod_surfcon.o mod_date.o
ini_stepone.o          : ini_stepone.f90 mod_tsteps.o
ini_fordate.o          : ini_fordate.f90 mod_dyncon0.o mod_lflags.o mod_atparam.o mod_hdifcon.o mod_physcon.o\
						 mod_surfcon.o mod_cli_land.o mod_date.o mod_var_land.o mod_cli_sea.o mod_var_sea.o\
						 mod_radcon.o
#land2atm.o             : land2atm.f
#land_init.o            : land_init.f
mod_atparam.o          : mod_atparam.f90
mod_cli_land.o         : mod_cli_land.f90 mod_atparam.o
mod_cli_sea.o          : mod_cli_sea.f90 mod_atparam.o
mod_cnvcon.o           : mod_cnvcon.f90
mod_cplcon_sea.o       : mod_cplcon_sea.f90 mod_atparam.o
mod_cplvar_sea.o       : mod_cplvar_sea.f90 mod_atparam.o
mod_cpl_flags.o        : mod_cpl_flags.f90
mod_date.o             : mod_date.f90 mod_atparam.o mod_tsteps.o
mod_dyncon0.o          : mod_dyncon0.f90
mod_dyncon1.o          : mod_dyncon1.f90 mod_atparam.o
mod_dyncon2.o          : mod_dyncon2.f90 mod_atparam.o
mod_dynvar.o           : mod_dynvar.f90 mod_atparam.o
mod_fft.o              : mod_fft.f90 mod_atparam.o
mod_flx_land.o         : mod_flx_land.f90 mod_atparam.o
mod_flx_sea.o          : mod_flx_sea.f90 mod_atparam.o
mod_hdifcon.o          : mod_hdifcon.f90 mod_atparam.o
mod_lflags.o           : mod_lflags.f90
mod_lsccon.o           : mod_lsccon.f90
mod_physcon.o          : mod_physcon.f90 mod_atparam.o
mod_physvar.o          : mod_physvar.f90 mod_atparam.o
mod_radcon.o           : mod_radcon.f90 mod_atparam.o
mod_randfor.o          : mod_randfor.f90 mod_atparam.o
mod_sflcon.o           : mod_sflcon.f90 mod_atparam.o
mod_spectral.o         : mod_spectral.f90 mod_atparam.o
mod_sppt.o             : mod_atparam.o mod_tsteps.o mod_dyncon1.o mod_spectral.o
mod_surfcon.o          : mod_surfcon.f90 mod_atparam.o
mod_tmean.o            : mod_tmean.f90 mod_atparam.o
mod_tsteps.o           : mod_tsteps.f90
mod_var_land.o         : mod_var_land.f90 mod_atparam.o
mod_var_sea.o          : mod_var_sea.f90 mod_atparam.o
mod_vdicon.o           : mod_vdicon.f90
phy_convmf.o           : phy_convmf.f90 mod_cnvcon.o mod_atparam.o mod_physcon.o
phy_lscond.o           : phy_lscond.f90 mod_lsccon.o mod_atparam.o mod_physcon.o
phy_phypar.o           : phy_phypar.f90 mod_cpl_flags.o mod_lflags.o mod_atparam.o mod_physcon.o mod_randfor.o \
					     mod_surfcon.o mod_var_land.o mod_cli_sea.o mod_var_sea.o mod_physvar.o mod_sppt.o\
                                             mod_tsteps.o
phy_radiat.o           : phy_radiat.f90 mod_atparam.o mod_physcon.o mod_radcon.o
phy_shtorh.o           : phy_shtorh.f90
phy_suflux.o           : phy_suflux.f90 mod_atparam.o mod_sflcon.o mod_physcon.o mod_radcon.o
phy_vdifsc.o           : phy_vdifsc.f90 mod_atparam.o mod_vdicon.o mod_physcon.o
ppo_diagns.o           : ppo_diagns.f90 mod_tsteps.o mod_atparam.o mod_dynvar.o mod_date.o
ppo_restart.o          : ppo_restart.f90 mod_tsteps.o mod_atparam.o mod_dynvar.o mod_date.o
ppo_setctl.o           : ppo_setctl.f90
ppo_setgrd.o           : ppo_setgrd.f90 mod_tsteps.o mod_date.o
ppo_tminc.o            : ppo_tminc.f90 mod_lflags.o mod_atparam.o mod_tmean.o mod_physcon.o mod_surfcon.o\
						 mod_cli_land.o mod_cli_sea.o mod_var_sea.o mod_physvar.o mod_radcon.o
ppo_dmflux.o           : ppo_dmflux.f90 mod_tsteps.o mod_atparam.o mod_tmean.o mod_flx_land.o mod_physcon.o\
						 mod_surfcon.o mod_flx_sea.o mod_var_sea.o mod_physvar.o mod_radcon.o
ppo_tmout.o            : ppo_tmout.f90 mod_tsteps.o mod_atparam.o mod_tmean.o mod_physcon.o
ppo_dmout.o            : ppo_dmout.f90 mod_tsteps.o mod_atparam.o mod_tmean.o
ppo_iogrid.o           : ppo_iogrid.f90 speedy_res_interface.o mod_io.o
#sea2atm.o              : sea2atm.f
#sea_init.o             : sea_init.f 
#send_to_coupler.o      : send_to_coupler.f
spe_matinv.o           : spe_matinv.f90
spe_spectral.o         : spe_spectral.f90 mod_atparam.o mod_spectral.o
spe_subfft_fftpack.o   : spe_subfft_fftpack.f90 mod_atparam.o mod_spectral.o mod_fft.o
