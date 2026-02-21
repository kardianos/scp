CC = gcc
CFLAGS = -O2 -Wall -Wextra -std=c11 -lm
LDFLAGS = -lm

BINDIR = ../bin
DATADIR = ../data

all: $(BINDIR)/sourced_maxwell $(BINDIR)/hessian_check $(BINDIR)/nullrotor_metric $(BINDIR)/depletion_gravity $(BINDIR)/monopole_coupling

$(BINDIR)/sourced_maxwell: sourced_maxwell.c clifford.h field.h
	@mkdir -p $(BINDIR)
	$(CC) $(CFLAGS) -o $@ sourced_maxwell.c $(LDFLAGS)

$(BINDIR)/hessian_check: hessian_check.c field.c clifford.h field.h
	@mkdir -p $(BINDIR)
	$(CC) $(CFLAGS) -fopenmp -o $@ hessian_check.c field.c $(LDFLAGS)

$(BINDIR)/nullrotor_metric: nullrotor_metric.c
	@mkdir -p $(BINDIR)
	$(CC) $(CFLAGS) -o $@ nullrotor_metric.c $(LDFLAGS)

$(BINDIR)/depletion_gravity: depletion_gravity.c field.c clifford.h field.h
	@mkdir -p $(BINDIR)
	$(CC) $(CFLAGS) -fopenmp -o $@ depletion_gravity.c field.c $(LDFLAGS)

$(BINDIR)/monopole_coupling: monopole_coupling.c
	@mkdir -p $(BINDIR)
	$(CC) $(CFLAGS) -o $@ monopole_coupling.c $(LDFLAGS)

run_monopole: $(BINDIR)/monopole_coupling
	@mkdir -p $(DATADIR)
	cd .. && bin/monopole_coupling -profile data/profile_sigma_e1.dat -outdir data

run_depletion: $(BINDIR)/depletion_gravity
	@mkdir -p $(DATADIR)
	cd .. && bin/depletion_gravity -profile data/profile_sigma_e1.dat -outdir data -mode A
	cd .. && bin/depletion_gravity -profile data/profile_sigma_e1.dat -outdir data -mode B
	cd .. && bin/depletion_gravity -profile data/profile_sigma_e1.dat -outdir data -mode C
	cd .. && bin/depletion_gravity -profile data/profile_sigma_e1.dat -outdir data -mode D

run_metric: $(BINDIR)/nullrotor_metric
	@mkdir -p $(DATADIR)
	cd .. && bin/nullrotor_metric -profile data/profile_sigma_e1.dat -outdir data

run_maxwell: $(BINDIR)/sourced_maxwell
	@mkdir -p $(DATADIR)
	cd .. && bin/sourced_maxwell -profile data/profile_sigma_e1.dat -outdir data

run_hessian: $(BINDIR)/hessian_check
	@mkdir -p $(DATADIR)
	cd .. && bin/hessian_check -profile data/profile_sigma_e1.dat -outdir data

clean:
	rm -f $(BINDIR)/sourced_maxwell $(BINDIR)/hessian_check $(BINDIR)/nullrotor_metric $(BINDIR)/depletion_gravity $(BINDIR)/monopole_coupling

.PHONY: all clean run_maxwell run_hessian run_metric run_depletion run_monopole
