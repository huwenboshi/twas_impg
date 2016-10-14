CC = g++
COPTS = -O3 -Wall

BINS = ImpG-Summary-GenBeta ImpG-Summary ImpG-Summary-Cor \
		ImpG-Summary-GenBeta-Cor ImpG-Summary-GenBeta-Tmp
all: $(BINS)

# ImpG-Summary
ImpG-Summary: impg.o impz_main.o util.o linsubs.o
	$(CC) impg.o impz_main.o util.o linsubs.o \
		-o ImpG-Summary $(COPTS)
impg_summary_main.o: impz_main.cpp
	$(CC) -c impz_main.cpp $(COPTS)

# ImpG-Summary-Cor
ImpG-Summary-Cor: impg.o impz_cor_main.o util.o linsubs.o
	$(CC) impg.o impz_cor_main.o util.o linsubs.o \
		-o ImpG-Summary-Cor $(COPTS)
impg_summary_cor_main.o: impz_cor_main.cpp
	$(CC) -c impz_cor_main.cpp $(COPTS)

# ImpG-Summary-GenBeta-Tmp
ImpG-Summary-GenBeta-Tmp: impg.o zgenbt_main_tmp.o util.o linsubs.o
	$(CC) impg.o zgenbt_main_tmp.o util.o linsubs.o \
		-o ImpG-Summary-GenBeta-Tmp $(COPTS)
zgenbt_main_tmp.o: zgenbt_main_tmp.cpp
	$(CC) -c zgenbt_main_tmp.cpp $(COPTS)

# ImpG-Summary-GenBeta
ImpG-Summary-GenBeta: impg.o zgenbt_main.o util.o linsubs.o
	$(CC) impg.o zgenbt_main.o util.o linsubs.o \
		-o ImpG-Summary-GenBeta $(COPTS)
zgenbt_main.o: zgenbt_main.cpp
	$(CC) -c zgenbt_main.cpp $(COPTS)

# ImpG-Summary-GenBeta-Cor
ImpG-Summary-GenBeta-Cor: impg.o zgenbt_cor_main.o util.o linsubs.o
	$(CC) impg.o zgenbt_cor_main.o util.o linsubs.o \
		-o ImpG-Summary-GenBeta-Cor $(COPTS)
zgenbt_cor_main.o: zgenbt_cor_main.cpp
	$(CC) -c zgenbt_cor_main.cpp $(COPTS)

# ImpG-Summary commons
impg.o: impg.cpp impg.h
	$(CC) -c impg.cpp $(COPTS)

# all common
util.o: util.cpp util.h
	$(CC) -c util.cpp $(COPTS)
linsubs.o: linsubs.cpp linsubs.h
	$(CC) -c linsubs.cpp $(COPTS)

clean:
	rm -f *.o $(BINS)
