CC = gcc
CFLAGS = -O3 -march=native -fopenmp -Wall -Wextra -std=c11
LDFLAGS = -lm -fopenmp
SRCDIR = src
BUILDDIR = build
BINDIR = bin

# Main targets
EVOLVE_OBJS = $(BUILDDIR)/hopfion_evolve.o $(BUILDDIR)/spherical_grid.o \
              $(BUILDDIR)/hopfion_field.o $(BUILDDIR)/hopfion_init.o \
              $(BUILDDIR)/topology.o
EVOLVE_TARGET = $(BINDIR)/hopfion_evolve

TEST_OBJS = $(BUILDDIR)/hopfion_test.o $(BUILDDIR)/spherical_grid.o \
            $(BUILDDIR)/hopfion_field.o $(BUILDDIR)/hopfion_init.o \
            $(BUILDDIR)/topology.o
TEST_TARGET = $(BINDIR)/hopfion_test

VERIFY_OBJS = $(BUILDDIR)/verify_gradient.o $(BUILDDIR)/spherical_grid.o \
              $(BUILDDIR)/hopfion_field.o $(BUILDDIR)/hopfion_init.o \
              $(BUILDDIR)/topology.o
VERIFY_TARGET = $(BINDIR)/verify_gradient

EQUIL_OBJS = $(BUILDDIR)/skyrmion_equil.o $(BUILDDIR)/spherical_grid.o \
             $(BUILDDIR)/hopfion_field.o $(BUILDDIR)/hopfion_init.o
EQUIL_TARGET = $(BINDIR)/skyrmion_equil

METRIC_OBJS = $(BUILDDIR)/hopfion_metric.o $(BUILDDIR)/spherical_grid.o \
              $(BUILDDIR)/hopfion_field.o $(BUILDDIR)/hopfion_init.o \
              $(BUILDDIR)/topology.o
METRIC_TARGET = $(BINDIR)/hopfion_metric

.PHONY: all clean evolve test verify equil metric

all: $(EVOLVE_TARGET)

evolve: $(EVOLVE_TARGET)

test: $(TEST_TARGET)
	./$(TEST_TARGET)

verify: $(VERIFY_TARGET)
	./$(VERIFY_TARGET)

equil: $(EQUIL_TARGET)
	./$(EQUIL_TARGET)

metric: $(METRIC_TARGET)
	./$(METRIC_TARGET)

$(EVOLVE_TARGET): $(EVOLVE_OBJS) | $(BINDIR)
	$(CC) $(EVOLVE_OBJS) -o $@ $(LDFLAGS)

$(TEST_TARGET): $(TEST_OBJS) | $(BINDIR)
	$(CC) $(TEST_OBJS) -o $@ $(LDFLAGS)

$(VERIFY_TARGET): $(VERIFY_OBJS) | $(BINDIR)
	$(CC) $(VERIFY_OBJS) -o $@ $(LDFLAGS)

$(EQUIL_TARGET): $(EQUIL_OBJS) | $(BINDIR)
	$(CC) $(EQUIL_OBJS) -o $@ $(LDFLAGS)

$(METRIC_TARGET): $(METRIC_OBJS) | $(BINDIR)
	$(CC) $(METRIC_OBJS) -o $@ $(LDFLAGS)

$(BUILDDIR)/%.o: $(SRCDIR)/%.c | $(BUILDDIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILDDIR):
	mkdir -p $(BUILDDIR)

$(BINDIR):
	mkdir -p $(BINDIR)

# Dependencies
$(BUILDDIR)/spherical_grid.o: $(SRCDIR)/spherical_grid.c $(SRCDIR)/spherical_grid.h $(SRCDIR)/clifford.h
$(BUILDDIR)/hopfion_field.o: $(SRCDIR)/hopfion_field.c $(SRCDIR)/hopfion_field.h $(SRCDIR)/spherical_grid.h $(SRCDIR)/clifford.h
$(BUILDDIR)/hopfion_init.o: $(SRCDIR)/hopfion_init.c $(SRCDIR)/hopfion_init.h $(SRCDIR)/spherical_grid.h $(SRCDIR)/clifford.h
$(BUILDDIR)/topology.o: $(SRCDIR)/topology.c $(SRCDIR)/topology.h $(SRCDIR)/hopfion_field.h $(SRCDIR)/spherical_grid.h $(SRCDIR)/clifford.h
$(BUILDDIR)/hopfion_evolve.o: $(SRCDIR)/hopfion_evolve.c $(SRCDIR)/spherical_grid.h $(SRCDIR)/hopfion_field.h $(SRCDIR)/hopfion_init.h $(SRCDIR)/topology.h $(SRCDIR)/clifford.h
$(BUILDDIR)/skyrmion_equil.o: $(SRCDIR)/skyrmion_equil.c $(SRCDIR)/spherical_grid.h $(SRCDIR)/hopfion_field.h $(SRCDIR)/hopfion_init.h $(SRCDIR)/clifford.h
$(BUILDDIR)/verify_gradient.o: $(SRCDIR)/verify_gradient.c $(SRCDIR)/spherical_grid.h $(SRCDIR)/hopfion_field.h $(SRCDIR)/hopfion_init.h $(SRCDIR)/clifford.h
$(BUILDDIR)/hopfion_metric.o: $(SRCDIR)/hopfion_metric.c $(SRCDIR)/spherical_grid.h $(SRCDIR)/hopfion_field.h $(SRCDIR)/hopfion_init.h $(SRCDIR)/topology.h $(SRCDIR)/clifford.h
$(BUILDDIR)/hopfion_test.o: $(SRCDIR)/hopfion_test.c $(SRCDIR)/spherical_grid.h $(SRCDIR)/hopfion_field.h $(SRCDIR)/hopfion_init.h $(SRCDIR)/topology.h $(SRCDIR)/clifford.h

clean:
	rm -rf $(BUILDDIR) $(BINDIR)
