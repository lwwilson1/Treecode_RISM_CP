OMPOFILES = treedriver.o tree_cp.o
OMPTARGET = treecp-omp

OFILES = treedriver.o tree_cp.o
TARGET = treecp

AOFILES = exact_coul.o 
ATARGET = direct

BOFILES = exact_coul_grid.o 
BTARGET = direct-grid

CC = mpifort
FFLAGS = -Wall -Wextra -O3 -qopenmp

$(OMPTARGET): $(OMPOFILES)
	$(CC) $(OMPOFILES) -qopenmp -o $@

$(TARGET): $(OFILES)
	$(CC) $(OFILES) -o $@

$(ATARGET): $(AOFILES)
	$(CC) $(AOFILES) -o $@

$(BTARGET): $(BOFILES)
	$(CC) $(BOFILES) -o $@


treedriver.o: treedriver.f90
	$(CC) $(FFLAGS) -c $<

tree_cp.o: tree_cp.f90 
	$(CC) $(FFLAGS) -c $<

exact_coul.o: exact_coul.f90
	$(CC) $(FFLAGS) -c $<

exact_coul_grid.o: exact_coul_grid.f90
	$(CC) $(FFLAGS) -c $<

treecode_procedures.mod: tree_cp.f90 tree_cp.o
	@true


clean:
	rm -f $(OFILES) $(AOFILES) $(BOFILES) $(TARGET) $(ATARGET) $(BTARGET) *.mod

