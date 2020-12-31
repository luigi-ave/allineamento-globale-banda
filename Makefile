main.exe: main.o needleman_wunsch_funct.o
	g++ main.o needleman_wunsch_funct.o -o main.exe

main.o: main.c
	g++ -c main.c -o main.o

needleman_wunsch_funct.o: needleman_wunsch_funct.c
	g++ -c needleman_wunsch_funct.c -o needleman_wunsch_funct.o

.PHONY: clean
clean:
	rm -r *.o *.exe