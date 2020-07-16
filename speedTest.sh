# start this script with command:
#   unbuffer ./speedTest.sh | tee -a ./speedGccBenchmark_CPU_SMT_OFF_j16.txt
# Note: unbuffer is in package 'expect'



## source ~/icc/bin/iccvars.sh intel64
#git lo -n 30
#git diff

date
#make uninstall && nice -n 40 make all JOBSNUM=16 DECI=6   MPFR=0 SSE=0 NAT=0 FETS=-8 THE_FAST=EXPERIMENTALLY_FAST=OFF    && echo "\e[96m==gcc  =A 1 ===\e[0m"                                                                                                  ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
#make uninstall && nice -n 40 make all JOBSNUM=16 DECI=6   MPFR=0 SSE=0 NAT=0 FETS=-8 THE_FAST=DANGEROUSLY_FAST=OFF SSE=1 && echo "\e[96m==gcc  =A 2 ===\e[0m"                                                                                                  ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
 make uninstall && nice -n 40 make all JOBSNUM=16 DECI=6   MPFR=0 SSE=0 NAT=0 FETS=-8 THE_FAST=EXPERIMENTALLY_FAST=ON     && echo "\e[96m==gcc  =A 3 ===\e[0m"                                                                                                  ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
#make uninstall && nice -n 40 make all JOBSNUM=16 DECI=6   MPFR=0 SSE=0 NAT=0 FETS=-8 THE_FAST=DANGEROUSLY_FAST=ON        && echo "\e[96m==gcc  =A 4 ===\e[0m"                                                                                                  ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
#make uninstall && nice -n 40 make all JOBSNUM=16 DECI=6   MPFR=0 SSE=0 NAT=0 FETS=-8 THE_FAST=DANGEROUSLY_FAST=ON  SSE=1 && echo "\e[96m==gcc  =A 5 ===\e[0m"                                                                                                  ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
date
#make uninstall && nice -n 40 make all JOBSNUM=16 DECI=15  MPFR=0 SSE=0 NAT=0 FETS=-8 THE_FAST=EXPERIMENTALLY_FAST=OFF    && echo "\e[96m==gcc  =B 1 ===\e[0m"                                                                                                  ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
#make uninstall && nice -n 40 make all JOBSNUM=16 DECI=15  MPFR=0 SSE=0 NAT=0 FETS=-8 THE_FAST=DANGEROUSLY_FAST=OFF SSE=1 && echo "\e[96m==gcc  =B 2 ===\e[0m"                                                                                                  ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
 make uninstall && nice -n 40 make all JOBSNUM=16 DECI=15  MPFR=0 SSE=0 NAT=0 FETS=-8 THE_FAST=EXPERIMENTALLY_FAST=ON     && echo "\e[96m==gcc  =B 3 ===\e[0m"                                                                                                  ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
#make uninstall && nice -n 40 make all JOBSNUM=16 DECI=15  MPFR=0 SSE=0 NAT=0 FETS=-8 THE_FAST=DANGEROUSLY_FAST=ON        && echo "\e[96m==gcc  =B 4 ===\e[0m"                                                                                                  ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
#make uninstall && nice -n 40 make all JOBSNUM=16 DECI=15  MPFR=0 SSE=0 NAT=0 FETS=-8 THE_FAST=DANGEROUSLY_FAST=ON  SSE=1 && echo "\e[96m==gcc  =B 5 ===\e[0m"                                                                                                  ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
date
#make uninstall && nice -n 40 make all JOBSNUM=16 DECI=18  MPFR=0 SSE=0 NAT=0 FETS=-8 THE_FAST=EXPERIMENTALLY_FAST=OFF    && echo "\e[96m==gcc  =C 1 ===\e[0m"                                                                                                  ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
#make uninstall && nice -n 40 make all JOBSNUM=16 DECI=18  MPFR=0 SSE=0 NAT=0 FETS=-8 THE_FAST=DANGEROUSLY_FAST=OFF SSE=1 && echo "\e[96m==gcc  =C 2 ===\e[0m"                                                                                                  ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
 make uninstall && nice -n 40 make all JOBSNUM=16 DECI=18  MPFR=0 SSE=0 NAT=0 FETS=-8 THE_FAST=EXPERIMENTALLY_FAST=ON     && echo "\e[96m==gcc  =C 3 ===\e[0m"                                                                                                  ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
#make uninstall && nice -n 40 make all JOBSNUM=16 DECI=18  MPFR=0 SSE=0 NAT=0 FETS=-8 THE_FAST=DANGEROUSLY_FAST=ON        && echo "\e[96m==gcc  =C 4 ===\e[0m"                                                                                                  ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
#make uninstall && nice -n 40 make all JOBSNUM=16 DECI=18  MPFR=0 SSE=0 NAT=0 FETS=-8 THE_FAST=DANGEROUSLY_FAST=ON  SSE=1 && echo "\e[96m==gcc  =C 5 ===\e[0m"                                                                                                  ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
date
#make uninstall && nice -n 40 make all JOBSNUM=16 DECI=33  MPFR=0 SSE=0 NAT=0 FETS=-8 THE_FAST=EXPERIMENTALLY_FAST=OFF    && echo "\e[96m==gcc  =D 1 ===\e[0m"                                                                                                  ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
#make uninstall && nice -n 40 make all JOBSNUM=16 DECI=33  MPFR=0 SSE=0 NAT=0 FETS=-8 THE_FAST=DANGEROUSLY_FAST=OFF SSE=1 && echo "\e[96m==gcc  =D 2 ===\e[0m"                                                                                                  ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
 make uninstall && nice -n 40 make all JOBSNUM=16 DECI=33  MPFR=0 SSE=0 NAT=0 FETS=-8 THE_FAST=EXPERIMENTALLY_FAST=ON     && echo "\e[96m==gcc  =D 3 ===\e[0m"                                                                                                  ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
#make uninstall && nice -n 40 make all JOBSNUM=16 DECI=33  MPFR=0 SSE=0 NAT=0 FETS=-8 THE_FAST=DANGEROUSLY_FAST=ON        && echo "\e[96m==gcc  =D 4 ===\e[0m"                                                                                                  ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
#make uninstall && nice -n 40 make all JOBSNUM=16 DECI=33  MPFR=0 SSE=0 NAT=0 FETS=-8 THE_FAST=DANGEROUSLY_FAST=ON  SSE=1 && echo "\e[96m==gcc  =D 5 ===\e[0m"                                                                                                  ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
date
#make uninstall && nice -n 40 make all JOBSNUM=16 DECI=62  MPFR=1 SSE=0 NAT=0 FETS=-8 THE_FAST=EXPERIMENTALLY_FAST=OFF    && echo "\e[96m==gcc  =E 1 ===\e[0m"                                                                                                  ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
#make uninstall && nice -n 40 make all JOBSNUM=16 DECI=62  MPFR=1 SSE=0 NAT=0 FETS=-8 THE_FAST=DANGEROUSLY_FAST=OFF SSE=1 && echo "\e[96m==gcc  =E 2 ===\e[0m"                                                                                                  ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
 make uninstall && nice -n 40 make all JOBSNUM=16 DECI=62  MPFR=1 SSE=0 NAT=0 FETS=-8 THE_FAST=EXPERIMENTALLY_FAST=ON     && echo "\e[96m==gcc  =E 3 ===\e[0m"                                                                                                  ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
#make uninstall && nice -n 40 make all JOBSNUM=16 DECI=62  MPFR=1 SSE=0 NAT=0 FETS=-8 THE_FAST=DANGEROUSLY_FAST=ON        && echo "\e[96m==gcc  =E 4 ===\e[0m"                                                                                                  ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
#make uninstall && nice -n 40 make all JOBSNUM=16 DECI=62  MPFR=1 SSE=0 NAT=0 FETS=-8 THE_FAST=DANGEROUSLY_FAST=ON  SSE=1 && echo "\e[96m==gcc  =E 5 ===\e[0m"                                                                                                  ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
date
#make uninstall && nice -n 40 make all JOBSNUM=16 DECI=150 MPFR=1 SSE=0 NAT=0 FETS=-8 THE_FAST=EXPERIMENTALLY_FAST=OFF    && echo "\e[96m==gcc  =F 1 ===\e[0m"                                                                                                  ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
#make uninstall && nice -n 40 make all JOBSNUM=16 DECI=150 MPFR=1 SSE=0 NAT=0 FETS=-8 THE_FAST=DANGEROUSLY_FAST=OFF SSE=1 && echo "\e[96m==gcc  =F 2 ===\e[0m"                                                                                                  ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
 make uninstall && nice -n 40 make all JOBSNUM=16 DECI=150 MPFR=1 SSE=0 NAT=0 FETS=-8 THE_FAST=EXPERIMENTALLY_FAST=ON     && echo "\e[96m==gcc  =F 3 ===\e[0m"                                                                                                  ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
#make uninstall && nice -n 40 make all JOBSNUM=16 DECI=150 MPFR=1 SSE=0 NAT=0 FETS=-8 THE_FAST=DANGEROUSLY_FAST=ON        && echo "\e[96m==gcc  =F 4 ===\e[0m"                                                                                                  ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
#make uninstall && nice -n 40 make all JOBSNUM=16 DECI=150 MPFR=1 SSE=0 NAT=0 FETS=-8 THE_FAST=DANGEROUSLY_FAST=ON  SSE=1 && echo "\e[96m==gcc  =F 5 ===\e[0m"                                                                                                  ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
date
#make uninstall && nice -n 40 make all JOBSNUM=16 DECI=62  MPFR=0 SSE=0 NAT=0 FETS=-8 THE_FAST=EXPERIMENTALLY_FAST=OFF    && echo "\e[96m==gcc  =G 1 ===\e[0m"                                                                                                  ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
#make uninstall && nice -n 40 make all JOBSNUM=16 DECI=62  MPFR=0 SSE=0 NAT=0 FETS=-8 THE_FAST=DANGEROUSLY_FAST=OFF SSE=1 && echo "\e[96m==gcc  =G 2 ===\e[0m"                                                                                                  ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
 make uninstall && nice -n 40 make all JOBSNUM=16 DECI=62  MPFR=0 SSE=0 NAT=0 FETS=-8 THE_FAST=EXPERIMENTALLY_FAST=ON     && echo "\e[96m==gcc  =G 3 ===\e[0m"                                                                                                  ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
#make uninstall && nice -n 40 make all JOBSNUM=16 DECI=62  MPFR=0 SSE=0 NAT=0 FETS=-8 THE_FAST=DANGEROUSLY_FAST=ON        && echo "\e[96m==gcc  =G 4 ===\e[0m"                                                                                                  ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
#make uninstall && nice -n 40 make all JOBSNUM=16 DECI=62  MPFR=0 SSE=0 NAT=0 FETS=-8 THE_FAST=DANGEROUSLY_FAST=ON  SSE=1 && echo "\e[96m==gcc  =G 5 ===\e[0m"                                                                                                  ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
date

 make uninstall && nice -n 40 make all JOBSNUM=16 DECI=6   MPFR=0 SSE=0 NAT=1 FETS=-8 THE_FAST=EXPERIMENTALLY_FAST=ON     && echo "\e[96m==gcc  =A 3 native===\e[0m"                                                                                            ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
 make uninstall && nice -n 40 make all JOBSNUM=16 DECI=15  MPFR=0 SSE=0 NAT=1 FETS=-8 THE_FAST=EXPERIMENTALLY_FAST=ON     && echo "\e[96m==gcc  =B 3 native===\e[0m"                                                                                            ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
 make uninstall && nice -n 40 make all JOBSNUM=16 DECI=18  MPFR=0 SSE=0 NAT=1 FETS=-8 THE_FAST=EXPERIMENTALLY_FAST=ON     && echo "\e[96m==gcc  =C 3 native===\e[0m"                                                                                            ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
 make uninstall && nice -n 40 make all JOBSNUM=16 DECI=33  MPFR=0 SSE=0 NAT=1 FETS=-8 THE_FAST=EXPERIMENTALLY_FAST=ON     && echo "\e[96m==gcc  =D 3 native===\e[0m"                                                                                            ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
 make uninstall && nice -n 40 make all JOBSNUM=16 DECI=62  MPFR=1 SSE=0 NAT=1 FETS=-8 THE_FAST=EXPERIMENTALLY_FAST=ON     && echo "\e[96m==gcc  =E 3 native===\e[0m"                                                                                            ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
 make uninstall && nice -n 40 make all JOBSNUM=16 DECI=150 MPFR=1 SSE=0 NAT=1 FETS=-8 THE_FAST=EXPERIMENTALLY_FAST=ON     && echo "\e[96m==gcc  =F 3 native===\e[0m"                                                                                            ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
 make uninstall && nice -n 40 make all JOBSNUM=16 DECI=62  MPFR=0 SSE=0 NAT=1 FETS=-8 THE_FAST=EXPERIMENTALLY_FAST=ON     && echo "\e[96m==gcc  =G 3 native===\e[0m"                                                                                            ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 

