gcc -Wall -Wno-unused-function -O2 -DHAVE_PTHREAD  -D_KSW_MAIN -o ksw ksw.c malloc_wrap.o  utils.o -lz
