d = Debug
FC = gfortran-8

LIBTEC = /home/abc/bin/lib
OBJ_PATH = .

#PGOPT = -pg

exe = ML

# список модулей - транслируются заранее:
srcs.mods := $(wildcard ./code/modules/*.*)
# список файлов .f90 из директории code:
srcs.f90 := $(wildcard ./code/*.f90)
# список файлов .f90 из директории code:
srcs.for := $(wildcard ./code/*.for)

# замена и директории, и расширения в списке srcs.f90:
objs.f90 := $(patsubst ./code/%.f90, %.o, $(srcs.f90))
# замена и директории, и расширения в списке srcs.for:
objs.for := $(patsubst ./code/%.for, %.o, $(srcs.for))
# замена и директории, и расширения в списке srcs.mods:
objs.mods := $(patsubst ./code/modules/%.f90, %.o, $(srcs.mods))

##mods := variables.for
##srcs := $(mods) main.for side.for cell.for init.for inout-sX.F90 inout-sY.F90 inout-s3D.F90 inout.for
#objs := main.o side.o cell.o
##objs := $(patsubst %.for, %.o,$(srcs))
##objs := $(patsubst %.F90, %.o,$(objs))

FCOPTS = -cpp -std=gnu -fdec-structure -fimplicit-none -finit-real=nan -fcray-pointer -ffree-line-length-0 -ffixed-line-length-0 -DR8=8

FCD = -g3 -O0 -DDLEV=1 -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow -Wno-array-temporaries
#FCR = -Ofast -march=native -floop-parallelize-all $(PGOPT)
FCR = -O3 $(PGOPT)
 
#FCOPTS = -Wl,-rpath,${CONDA_PREFIX}/lib -I ~/bin/lib
#FCOPTS = -O3 -fdefault-real-8 -std=gnu -fdec-structure -fimplicit-none -finit-real=nan -fcray-pointer

ifeq ($(DEB),1)
  FFLAGS = $(FCOPTS) $(FCD)
else
  FFLAGS = $(FCOPTS) $(FCR)
endif

FLOPTS = $(PGOPT) -L ./lib -ltecio #-lpthread -lstdc++ -lc++
#LDFLAGS = -g

#$(info srcs.mods: ${srcs.mods})
#$(info srcs.f90: ${srcs.f90})
#$(info objs.mods: ${objs.mods})
#$(info FFLAGS: $(FFLAGS))

all: $(objs.mods) $(exe)

$(exe): $(objs.mods) $(objs.for) $(objs.f90)
	$(FC) $(PGOPT) $(LDFLAGS) $(objs.for) $(objs.f90) $(objs.mods) $(LIBTEC)/libtecio.a -lstdc++ -lpthread -o $(exe)

$(objs.mods): $(OBJ_PATH)/%.o : code/modules/%.f90
	$(FC) $(FFLAGS) -c $< -o $(OBJ_PATH)/$(@F) -I. -I./$(OBJ_PATH)/  -J $(OBJ_PATH)/

$(objs.f90): $(srcs.mods)
$(OBJ_PATH)/%.o: code/%.f90
	$(FC) $(FFLAGS) -c $< -o $(OBJ_PATH)/$(@F) -I. -I./$(OBJ_PATH)/  -J $(OBJ_PATH)/

$(objs.for): $(srcs.mods)
$(OBJ_PATH)/%.o: code/%.for
	$(FC) $(FFLAGS) -c $< -o $(OBJ_PATH)/$(@F) -I. -I./$(OBJ_PATH)/  -J $(OBJ_PATH)/

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

echo:
	@echo "exe="$(exe)
	@echo "srcs.f90="$(srcs.f90)
	@echo "srcs.mods="$(srcs.mods)
	@echo "objs.mods="$(objs.mods)
	@echo "objs.for="$(objs.for)
	@echo "objs.f90="$(objs.f90)

clean:
	rm -f ./*.o
	rm -f ./*.mod
	rm -f ./*.exe
