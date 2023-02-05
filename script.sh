for reynolds in 500. 750. 1300. 1600.
do
	mkdir re_$reynolds 
	cp -r NS2D_temp.c streamdiff.h re_$reynolds
	cd re_$reynolds 
	sed s/nb_reynolds/$reynolds/g NS2D_temp.c > NS2D.c
	qcc -O2 -Wall -o NS2D NS2D.c -lm 
	./NS2D 
	cd .. 
done 	
