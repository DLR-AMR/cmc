if ENABLE_FORTRAN
# Clean up modules files in the root directory (if no module directory has been specified)
CLEANFILES += *.$(FC_MODEXT)

# Define a variable holding the module directory (for a rule below)
cmc_current_moddir = 

if WITH_MODDIR
# Updates for the modue output and include path (if a seperate module directory has been specified)
AM_FCFLAGS += $(FC_MODOUT)@CMC_MODDIR@ $(FC_MODINC)@CMC_MODDIR@
AM_CPPFLAGS += -I@CMC_MODDIR@

# Clean the module files in this directory
CLEANFILES += @CMC_MODDIR@/*.$(FC_MODEXT)

# Add the creation of the module directory as an order only prerequisite to the Fortran files
$(MODSOURCES): %.f90 : | create-moddir

.PHONY : create-moddir

# Rule to create the module directory
create-moddir:
	$(MKDIR_P) @CMC_MODDIR@

# Save the module directory
cmc_current_moddir += @CMC_MODDIR@/

endif

# If the install target is made, we will copy the module files into the include directory (after the installation of the header files)
install-data-hook:
	@cp -fp $(cmc_current_moddir)*.$(FC_MODEXT) $(includedir)

# Define dependencies of the Fortran modules (in case they depend on other modules)
# This needs to be done in order to ensure the correct build process in any case
# ...

# Define dependencies for all Fortran programs of the Fortran modules
# This needs to be done in order to ensure the correct build process in any case
# The example 'cmc_fortran_lossy_amr_compression' depends on the modules: cmc_fortron, cmc_netcdf_interface, cmc_amr_compression_interface
example/lossy_amr_nc_compression/cmc_fortran_lossy_amr_compression.o : src/fortran_interfaces/cmc_fortran.o src/fortran_interfaces/cmc_netcdf_interface.o src/fortran_interfaces/cmc_amr_compression_interface.o

endif
