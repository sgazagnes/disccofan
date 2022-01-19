
#Compiler and Linker
CC          := mpicc
CCMISC      := gcc 
#The Target Binary Program
TARGETA      := disccofan
TARGETB      := disccofan_interactive

#The Directories, Source, Includes, Objects, Binary and Resources
SRCDIR      := src/dct
SRCMISC     := src/misc
INCDIR      := 
BUILDDIR    := obj
TARGETDIR   := ./
RESDIR      := res
SRCEXT      := c
DEPEXT      := d
OBJEXT      := o

#Flags, Libraries and Includes
CFLAGS      :=  -O3  -fopenmp -pedantic -Wall -Wextra -Wundef -Wshadow  -Wcast-align -Wstrict-prototypes -Wwrite-strings -Wcast-qual -Wswitch-default -Wunreachable-code -Wformat=2 -Winit-self -march=native -g -std=gnu99  -Wno-pointer-arith
LIB         := -fopenmp -lm -lfreeimage -lhdf5 -lcfitsio  -L/usr/local/hdf5/lib/
LIBM        :=  -lgsl -lgslcblas  -lfftw3f_omp -lfftw3f
INC         := -I/usr/local/include -I/usr/local/hdf5/include -I/usr/include/cfitsio/
INCDEP      := 

#---------------------------------------------------------------------------------
#DO NOT EDIT BELOW THIS LINE
#---------------------------------------------------------------------------------
SOURCESA     := $(shell find $(SRCDIR) -type f \( -name *.$(SRCEXT) ! -name "*step_main.c" ! -name "moschini*" \))
SOURCESB     := $(shell find $(SRCDIR) -type f \( -name *.$(SRCEXT) ! -name "seq_main.c" ! -name "moschini*" \))
OBJECTSA     := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCESA:.$(SRCEXT)=.$(OBJEXT)))
OBJECTSB     := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCESB:.$(SRCEXT)=.$(OBJEXT)))

#Defauilt Make
all: directories  $(TARGETA) $(TARGETB) 

misc: polyImage check  
#Remake
remake: cleaner all

#Copy Resources from Resources Directory to Target Directory
#resources: directories
#	@cp $(RESDIR)/* $(TARGETDIR)/

#Make the Directories
directories:
	@mkdir -p $(BUILDDIR)

#Clean only Objecst
clean:
	@$(RM) -rf $(BUILDDIR)

#Full Clean, Objects and Binaries
cleaner: clean
	@$(RM) $(TARGETA) $(TARGETB) polyImage check bubbles kde

#Pull in dependency info for *existing* .o files
-include $(OBJECTSA:.$(OBJEXT)=.$(DEPEXT))
-include $(OBJECTSB:.$(OBJEXT)=.$(DEPEXT))


#Link
$(TARGETA): $(OBJECTSA)
	$(CC) -o $(TARGETDIR)/$(TARGETA) $^ $(LIB)

$(TARGETB): $(OBJECTSB)
	$(CC) -o $(TARGETDIR)/$(TARGETB) $^ $(LIB)


#Compile
$(BUILDDIR)/%.$(OBJEXT): $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(dir $@)
#	$(CC) $(CFLAGS) $(INC) -c -o $@ $<
	@$(CC) $(CFLAGS)  $(INC) -c $< -o $@
	@$(CC) $(CFLAGS) $(INCDEP) -MM $(SRCDIR)/$*.$(SRCEXT) > $(BUILDDIR)/$*.$(DEPEXT)
	@cp -f $(BUILDDIR)/$*.$(DEPEXT) $(BUILDDIR)/$*.$(DEPEXT).tmp
	@sed -e 's|.*:|$(BUILDDIR)/$*.$(OBJEXT):|' < $(BUILDDIR)/$*.$(DEPEXT).tmp > $(BUILDDIR)/$*.$(DEPEXT)
	@sed -e 's/.*://' -e 's/\\$$//' < $(BUILDDIR)/$*.$(DEPEXT).tmp | fmt -1 | sed -e 's/^ *//' -e 's/$$/:/' >> $(BUILDDIR)/$*.$(DEPEXT)
	@rm -f $(BUILDDIR)/$*.$(DEPEXT).tmp

bubbles: $(SRCMISC)/bubbles_an.c  $(SRCMISC)/eispack.h  $(SRCMISC)/eispack.c $(SRCMISC)/common.h
	$(CCMISC)  $(INC) -o $@ $^  $(LIB)

polyImage: $(SRCMISC)/polyImage.c $(SRCMISC)/misc.h $(SRCMISC)/misc.c $(SRCMISC)/common.h $(SRCMISC)/kmeans.h  $(SRCMISC)/kmeans.c  $(SRCMISC)/io.h  $(SRCMISC)/io.c  $(SRCMISC)/eispack.h $(SRCMISC)/eispack.c
	$(CCMISC)  $(INC) -o $@ $^  $(LIB) $(LIBM)

check: $(SRCMISC)/check.c
	$(CCMISC) $(CFLAGS) $(INC) -o $@ $^  $(LIB)

kde: $(SRCMISC)/kerneldensity.c $(SRCMISC)/kerneldensity.h
	$(CCMISC) $(CFLAGS) $(INC) -o $@ $^  $(LIB)

#Non-File Targets
.PHONY: all remake clean cleaner resources
