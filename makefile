all:gsdfmain pickphase readlogfile fitphasev fitphasev2

gsdfmain:gsdfmain.c createdatabase.h datastruct.h stafunction.h pathgrading.h parameter.h
	gcc -m32 -o gsdfmain gsdfmain.c -g sacio.a 

pickphase:pickphase.c parameter.h
	gcc -m32 -o pickphase pickphase.c -g sacio.a

readlogfile:readlogfile.c parameter.h
	gcc -m32 -o readlogfile readlogfile.c -g

fitphasev:fitphasev.c
	gcc -o fitphasev fitphasev.c -g -lm -lgsl -lgslcblas

fitphasev2:fitphasev2.c
	gcc -o fitphasev2 fitphasev2.c -g -lm -lgsl -lgslcblas
