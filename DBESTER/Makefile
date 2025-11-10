FC		= mpifrtpx
CC		= mpifccpx

FFLAGS		= -Kfast -g -Ntl_trt -Nquickdbg
CFLAGS		= -Kfast -g -Ntl_trt -Nquickdbg
LIBDIR		= ./LIBEMSES
LDFLAGS		= -L/center/local/apl/fx/hdf5-FX100-mpi/lib -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lz \
		  -L/center/local/apl/fx/fftw-3.3.4-fx100 -lfftw3

HDF5DIR		= /center/local/apl/fx/hdf5-FX100-mpi/include
FFTWDIR		= /center/local/apl/fx/fftw-3.3.4-fx100/include

LINKER		= $(FC)

COMMONHDRS	= oh_dim.h oh_stats.h
OHHDRS		= ohhelp1.h ohhelp2.h ohhelp3.h oh_part.h
OHOBJS		= ohhelp1.o ohhelp2.o ohhelp3.o

FHDRS		= ohhelp_f.h
FMODS		= oh_type.o   oh_mod1.o   oh_mod2.o   oh_mod3.o

KMODS		= paramt.o interf.o allcom.o namels.o hdf.o ssl2ex.o hdf5io.o wrapper.o lemses.o

FOBJS		= main.o emses.o esses.o \
		  input.o ohinit.o inital.o inipcl.o inipin.o \
		  inifld.o getccl.o capmtx.o chgntr.o \
		  rescal.o chkprm.o \
		  digset.o wrthdf.o hdfdig.o \
		  psolve1.o psolve2.o positn.o psuper1.o psuper2.o \
		  charge.o rhoupd.o \
		  efield.o bfield.o esfild.o surchg.o ptnvar.o \
		  grdcrt.o bdyfld.o rflfld.o fsmask.o \
		  trnscs.o poimix.o poirft.o \
		  poisnt.o poicst.o phbond.o \
		  frmttd.o spsort.o snap.o nborps.o chkres.o

COBJS		= mt19937arm.o psort-struct.o

OBJS		= $(FOBJS) $(COBJS) $(FMODS) $(KMODS) $(LMODS) $(OHOBJS)

mpiemses3D:	$(OBJS)
		$(LINKER) $(OBJS) $(FFLAGS) $(LDFLAGS) -o $@

$(FOBJS):%.o:	%.F90 $(FMODS) $(KMODS) $(LMODS) $(FHDRS) $(COMMONHDRS)
		$(FC) $(FFLAGS) -I$(FFTWDIR) -I$(HDF5DIR) -c $< -o $@
$(COBJS):%.o:	%.c $(OHHDRS)
		$(CC) $(CFLAGS) -c $< -o $@
%simulator.o:	sample.o

$(FMODS):%.o:	%.F90 $(COMMONHDRS)
		$(FC) $(FFLAGS) -I$(FFTWDIR) -I$(HDF5DIR) -c $< -o $@
$(KMODS):%.o:	%.F90
		$(FC) $(FFLAGS) -I$(FFTWDIR) -I$(HDF5DIR) -c $< -o $@
oh_mod1.o:	oh_type.o
oh_mod2.o:	oh_mod1.o
oh_mod3.o:	oh_mod2.o
allcom.o:	oh_type.o paramt.o lemses.o
lemses.o:	hdf5io.o wrapper.o
wrapper.o:	ssl2ex.o
bdyfld.o:	vical.fnc
psuper1.o:	setpos.fnc
psuper2.o:	setpos.fnc
$(OHOBJS):%.o:	%.c $(COMMONHDRS) $(OHHDRS)
		$(CC) $(CFLAGS) -c $< -o $@

clean:;
		rm *.o
clear:;
		mv hdf5io.o hdf5io.tmp; rm *.o; mv hdf5io.tmp hdf5io.o
