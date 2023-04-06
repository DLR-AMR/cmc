if ENABLE_FORTRAN
# Clean up modules files in the root directory (if no module directory has been specified)
CLEANFILES += *.$(FC_MODEXT)

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

else
# If no module directory is given, we will copy the module files in the include directory (after the installation of the header files)
install-data-hook:
	@cp -fp *.$(FC_MODEXT) $(includedir)

endif

endif