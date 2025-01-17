all: fcc_S1.inp

check: fcc_S1.inp
	cp Dusch_S1.dat eldip  fcc_S1.inp Freq_S0.dat  Freq_S1.dat  GradQ_S1.dat FCC -v
	cd FCC; fcclasses3 fcc_S1.inp
	xmgrace Quantics/p4+/spec_LS_TD_hwhm0.010.dat FCC/spec_LS_TD.dat

fcc_S1.inp:
	../../Scripts/op2fcc.py pyr4+.op

clean:
	@rm Dusch_S1.dat  Dusch_S2.dat  eldip  fcc_S1.inp  fcc_S2.inp  Freq_S0.dat  Freq_S1.dat  Freq_S2.dat  GradQ_S1.dat  GradQ_S2.dat -v; sleep 0
	@rm FCC/* -v; sleep 0
