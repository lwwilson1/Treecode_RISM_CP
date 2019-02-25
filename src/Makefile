OFILES = treedriver.o tree_cp.o 
TARGET = treecp

AOFILES = exact_coul.o 
ATARGET = direct

AOFILES = exact_coul_grid.o 
ATARGET = direct-grid

CC = gfortran
FFLAGS = -Wall -Wextra -O3

$(TARGET): $(OFILES)
	$(CC) $(OFILES) -o $@

$(ATARGET): $(AOFILES)
	$(CC) $(AOFILES) -o $@


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
	rm -f $(OFILES) $(AOFILES) $(TARGET) $(ATARGET) *.mod

