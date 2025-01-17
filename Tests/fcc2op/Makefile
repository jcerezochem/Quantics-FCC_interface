all: oper.op

check: oper.op
	cd Quantics; rm p4+ -rf; sleep 0; quantics -mnd p4+_VH.inp
	cd Quantics/p4+; reconvolute_TD -f auto -E0 -1
	xmgrace Quantics/p4+/spec_LS_TD_hwhm0.010.dat FCC/spec_LS_TD.dat

oper.op:
	../../Scripts/fcc2op.py fcc_S1.out

clean:
	@rm oper.op -v; sleep 0
	@rm Quantics/p4+ -rfv; sleep 0
	@rm Quantics/fort.21 -v; sleep 0
