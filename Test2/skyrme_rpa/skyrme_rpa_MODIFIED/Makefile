FC = gfortran

LEGACY = -std=legacy

PREPROC = -cpp

MCMODEL = -mcmodel=medium

OPT = $(PREPROC) -O3 -march=native

DBG = $(PREPROC) -D DEBUG -Og -g -Wall                              \
                                 -Wno-tabs                          \
                                 -Wno-compare-reals                 \
                                 -Wno-unused-label                  \
                                 -Wno-unused-dummy-argument         \
                                 -Wno-unused-parameter              \
                                 -Wextra                            \
                                 -Wconversion                       \
                                 -Warray-temporaries                \
                                 -Wsurprising                       \
                                 -Wline-truncation                  \
                                 -Wcharacter-truncation             \
                                 -fcheck=all                        \
                                 -fbacktrace                        \
#                                -ffpe-trap=zero,overflow,underflow \
#                                -fimplicit-none                    \
#                                -Wpedantic                         \

SRC = ./skyrme_rpa.f

run: $(SRC)
	$(FC) $(MCMODEL) $(OPT) $(LEGACY) $(SRC) -o run

original: $(SRC)
	$(FC) $(MCMODEL) -DORIGINAL $(OPT) $(LEGACY) $(SRC) -o run

dbg: $(SRC)
	$(FC) $(MCMODEL) $(DBG) $(LEGACY) $(SRC) -o dbg

clear:
	rm -f run dbg ./*.dat ./*.out ./out_matvec/*.mat ./out_matvec/*.vec 
