#comp = gfortran
#flag = 
#-g -fcheck=all -ffpe-trap=zero,overflow,underflow \
       #-Warray-temporaries -fbounds-check -O2
		
 comp = ifort
# 
 flag = -check all -traceback
		
		
		
executavel = main





dependencias =  mod_rndgen_multiple.f90 geraRRN.f90 modelosEpidemicosOGA.f90 main.f90

#dependencias =  mod_rndgen_multiple.f90 mod_lpa_comp.f90 main_medias.f90
#executavel = main_medias

default:
	$(comp) $(flag) $(dependencias) -o $(executavel)
    
clean:
	rm -f *.mod
