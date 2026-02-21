CC = gcc
CFLAGS = -O3 -march=native -fopenmp -Wall -Wextra -std=c11
LDFLAGS = -lm -fopenmp
SRCDIR = src
BUILDDIR = build
BINDIR = bin

SRCS = $(SRCDIR)/main.c $(SRCDIR)/field.c $(SRCDIR)/initial.c
OBJS = $(BUILDDIR)/main.o $(BUILDDIR)/field.o $(BUILDDIR)/initial.o
TARGET = $(BINDIR)/soliton_search

VERIFY_OBJS = $(BUILDDIR)/verify.o $(BUILDDIR)/field.o $(BUILDDIR)/initial.o
VERIFY_TARGET = $(BINDIR)/verify_gradient

RADIAL_TARGET = $(BINDIR)/radial_solver
RATMAP_TARGET = $(BINDIR)/rational_map_solver
FINLAM_TARGET = $(BINDIR)/finite_lambda_solver
V3D_OBJS = $(BUILDDIR)/verify3d.o $(BUILDDIR)/field.o $(BUILDDIR)/initial.o
V3D_TARGET = $(BINDIR)/verify3d
MAXWELL_TARGET = $(BINDIR)/maxwell
VEFF_TARGET = $(BINDIR)/veff
SCATTER_OBJS = $(BUILDDIR)/scatter.o $(BUILDDIR)/field.o $(BUILDDIR)/initial.o $(BUILDDIR)/coupling.o
SCATTER_TARGET = $(BINDIR)/scatter
VCOUP_OBJS = $(BUILDDIR)/verify_coupling.o $(BUILDDIR)/field.o $(BUILDDIR)/initial.o $(BUILDDIR)/coupling.o
VCOUP_TARGET = $(BINDIR)/verify_coupling
BOUND_TARGET = $(BINDIR)/bound
NMODES_TARGET = $(BINDIR)/normal_modes
AMODES_TARGET = $(BINDIR)/angular_modes
DEGEN_TARGET = $(BINDIR)/degenerate
VARSIG_TARGET = $(BINDIR)/varsig
TENSOR_TARGET = $(BINDIR)/tensor_prop
ONELOOP_TARGET = $(BINDIR)/oneloop
FPRIM_TARGET = $(BINDIR)/field_primary
PROTON_TARGET = $(BINDIR)/proton
HOPFION_TARGET = $(BINDIR)/hopfion
TORSION_TARGET = $(BINDIR)/torsion

.PHONY: all clean verify radial ratmap finlam v3d maxwell veff scatter vcoup bound nmodes amodes degen varsig tensor oneloop fprim proton hopfion torsion

all: $(TARGET) $(VERIFY_TARGET) $(RADIAL_TARGET) $(RATMAP_TARGET) $(FINLAM_TARGET) $(V3D_TARGET) $(MAXWELL_TARGET) $(VEFF_TARGET) $(SCATTER_TARGET) $(VCOUP_TARGET) $(BOUND_TARGET) $(NMODES_TARGET) $(AMODES_TARGET) $(DEGEN_TARGET) $(VARSIG_TARGET) $(TENSOR_TARGET) $(FPRIM_TARGET)

verify: $(VERIFY_TARGET)
	./$(VERIFY_TARGET)

radial: $(RADIAL_TARGET)
	./$(RADIAL_TARGET)

ratmap: $(RATMAP_TARGET)
	./$(RATMAP_TARGET)

finlam: $(FINLAM_TARGET)
	./$(FINLAM_TARGET)

v3d: $(V3D_TARGET)
	./$(V3D_TARGET)

$(TARGET): $(OBJS) | $(BINDIR)
	$(CC) $(OBJS) -o $@ $(LDFLAGS)

$(VERIFY_TARGET): $(VERIFY_OBJS) | $(BINDIR)
	$(CC) $(VERIFY_OBJS) -o $@ $(LDFLAGS)

$(RADIAL_TARGET): $(BUILDDIR)/radial.o | $(BINDIR)
	$(CC) $(BUILDDIR)/radial.o -o $@ -lm

$(RATMAP_TARGET): $(BUILDDIR)/rational_map.o | $(BINDIR)
	$(CC) $(BUILDDIR)/rational_map.o -o $@ -lm

$(FINLAM_TARGET): $(BUILDDIR)/finite_lambda.o | $(BINDIR)
	$(CC) $(BUILDDIR)/finite_lambda.o -o $@ -lm

$(V3D_TARGET): $(V3D_OBJS) | $(BINDIR)
	$(CC) $(V3D_OBJS) -o $@ $(LDFLAGS)

$(MAXWELL_TARGET): $(BUILDDIR)/maxwell.o | $(BINDIR)
	$(CC) $(BUILDDIR)/maxwell.o -o $@ -lm

$(VEFF_TARGET): $(BUILDDIR)/veff.o | $(BINDIR)
	$(CC) $(BUILDDIR)/veff.o -o $@ -lm

$(SCATTER_TARGET): $(SCATTER_OBJS) | $(BINDIR)
	$(CC) $(SCATTER_OBJS) -o $@ $(LDFLAGS)

$(VCOUP_TARGET): $(VCOUP_OBJS) | $(BINDIR)
	$(CC) $(VCOUP_OBJS) -o $@ $(LDFLAGS)

vcoup: $(VCOUP_TARGET)
	./$(VCOUP_TARGET)

bound: $(BOUND_TARGET)
	./$(BOUND_TARGET)

$(BOUND_TARGET): $(BUILDDIR)/bound.o | $(BINDIR)
	$(CC) $(BUILDDIR)/bound.o -o $@ -lm

nmodes: $(NMODES_TARGET)
	./$(NMODES_TARGET)

$(NMODES_TARGET): $(BUILDDIR)/normal_modes.o | $(BINDIR)
	$(CC) $(BUILDDIR)/normal_modes.o -o $@ -lm

amodes: $(AMODES_TARGET)
	./$(AMODES_TARGET)

$(AMODES_TARGET): $(BUILDDIR)/angular_modes.o | $(BINDIR)
	$(CC) $(BUILDDIR)/angular_modes.o -o $@ -lm

degen: $(DEGEN_TARGET)
	./$(DEGEN_TARGET)

$(DEGEN_TARGET): $(BUILDDIR)/degenerate.o | $(BINDIR)
	$(CC) $(BUILDDIR)/degenerate.o -o $@ -lm

varsig: $(VARSIG_TARGET)
	./$(VARSIG_TARGET)

$(VARSIG_TARGET): $(BUILDDIR)/varsig.o | $(BINDIR)
	$(CC) $(BUILDDIR)/varsig.o -o $@ -lm

tensor: $(TENSOR_TARGET)
	./$(TENSOR_TARGET)

$(TENSOR_TARGET): $(BUILDDIR)/tensor_prop.o | $(BINDIR)
	$(CC) $(BUILDDIR)/tensor_prop.o -o $@ -lm

oneloop: $(ONELOOP_TARGET)
	./$(ONELOOP_TARGET) -profile data/profiles/profile_sigma_e1.dat

$(ONELOOP_TARGET): $(BUILDDIR)/oneloop.o | $(BINDIR)
	$(CC) $(BUILDDIR)/oneloop.o -o $@ -lm

fprim: $(FPRIM_TARGET)
	./$(FPRIM_TARGET) -profile data/profiles/profile_sigma_e1.dat

$(FPRIM_TARGET): $(BUILDDIR)/field_primary.o | $(BINDIR)
	$(CC) $(BUILDDIR)/field_primary.o -o $@ -lm -fopenmp

proton: $(PROTON_TARGET)
	./$(PROTON_TARGET) -mode hedgehog -profile data/profiles/profile_sigma_e1.dat

$(PROTON_TARGET): $(BUILDDIR)/proton.o | $(BINDIR)
	$(CC) $(BUILDDIR)/proton.o -o $@ -lm -fopenmp

hopfion: $(HOPFION_TARGET)
	./$(HOPFION_TARGET) -N 64 -L 16 -relax 500 -steps 0

$(HOPFION_TARGET): $(BUILDDIR)/hopfion.o | $(BINDIR)
	$(CC) $(BUILDDIR)/hopfion.o -o $@ -lm -fopenmp

torsion: $(TORSION_TARGET)
	./$(TORSION_TARGET)

$(TORSION_TARGET): $(BUILDDIR)/torsion.o | $(BINDIR)
	$(CC) $(BUILDDIR)/torsion.o -o $@ -lm

$(BUILDDIR)/%.o: $(SRCDIR)/%.c | $(BUILDDIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILDDIR):
	mkdir -p $(BUILDDIR)

$(BINDIR):
	mkdir -p $(BINDIR)

# Dependencies
$(BUILDDIR)/main.o: $(SRCDIR)/main.c $(SRCDIR)/field.h $(SRCDIR)/initial.h $(SRCDIR)/clifford.h
$(BUILDDIR)/field.o: $(SRCDIR)/field.c $(SRCDIR)/field.h $(SRCDIR)/clifford.h
$(BUILDDIR)/initial.o: $(SRCDIR)/initial.c $(SRCDIR)/initial.h $(SRCDIR)/field.h $(SRCDIR)/clifford.h
$(BUILDDIR)/verify.o: $(SRCDIR)/verify.c $(SRCDIR)/field.h $(SRCDIR)/initial.h $(SRCDIR)/clifford.h
$(BUILDDIR)/radial.o: $(SRCDIR)/radial.c
$(BUILDDIR)/rational_map.o: $(SRCDIR)/rational_map.c
$(BUILDDIR)/finite_lambda.o: $(SRCDIR)/finite_lambda.c
$(BUILDDIR)/verify3d.o: $(SRCDIR)/verify3d.c $(SRCDIR)/field.h $(SRCDIR)/initial.h $(SRCDIR)/clifford.h
$(BUILDDIR)/maxwell.o: $(SRCDIR)/maxwell.c
$(BUILDDIR)/bound.o: $(SRCDIR)/bound.c
$(BUILDDIR)/normal_modes.o: $(SRCDIR)/normal_modes.c
$(BUILDDIR)/angular_modes.o: $(SRCDIR)/angular_modes.c
$(BUILDDIR)/degenerate.o: $(SRCDIR)/degenerate.c
$(BUILDDIR)/varsig.o: $(SRCDIR)/varsig.c
$(BUILDDIR)/oneloop.o: $(SRCDIR)/oneloop.c
$(BUILDDIR)/tensor_prop.o: $(SRCDIR)/tensor_prop.c
$(BUILDDIR)/field_primary.o: $(SRCDIR)/field_primary.c
$(BUILDDIR)/proton.o: $(SRCDIR)/proton.c
$(BUILDDIR)/hopfion.o: $(SRCDIR)/hopfion.c
$(BUILDDIR)/torsion.o: $(SRCDIR)/torsion.c
$(BUILDDIR)/veff.o: $(SRCDIR)/veff.c
$(BUILDDIR)/coupling.o: $(SRCDIR)/coupling.c $(SRCDIR)/coupling.h $(SRCDIR)/field.h $(SRCDIR)/clifford.h
$(BUILDDIR)/scatter.o: $(SRCDIR)/scatter.c $(SRCDIR)/field.h $(SRCDIR)/coupling.h $(SRCDIR)/clifford.h
$(BUILDDIR)/verify_coupling.o: $(SRCDIR)/verify_coupling.c $(SRCDIR)/field.h $(SRCDIR)/coupling.h $(SRCDIR)/initial.h $(SRCDIR)/clifford.h

clean:
	rm -rf $(BUILDDIR) $(BINDIR)
