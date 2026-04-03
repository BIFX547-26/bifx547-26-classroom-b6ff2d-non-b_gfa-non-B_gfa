PGMNAME = gfa
CC = gcc
OBJ = src/$(PGMNAME).o \
src/globals.o src/gfa_api.o \
src/cdna.o     src/findIR.o    src/nulls.o           src/process_repeats.o \
src/findAPR.o  src/findMR.o    src/print_gff_file.o  src/rcdna.o \
src/findDR.o   src/findSTR.o   src/is_subset.o  src/print_tsv_file.o  src/read_fasta.o \
src/findGQ.o   src/findZDNA.o  src/print_usage.o     src/read_mult_fasta.o
CFLAGS = -O2 -Isrc
LFLAGS = -lm
$(PGMNAME): $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) -o $@ $(LFLAGS)

clean:
	rm -f src/*.o $(PGMNAME)
