default:
	gcc -pthread pma.c -o run
	./run
clear:
	rm -rf out.exe run