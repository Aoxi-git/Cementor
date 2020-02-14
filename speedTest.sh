# start this script with command:
#   unbuffer ./speedTest.sh | tee -a ./speedIntelBenchmark.txt
# Note: unbuffer is in package 'expect'

# import Intel compiler settings and PATH
source ~/icc/bin/iccvars.sh intel64

#git lo -n 30
#git diff

date
make uninstall && nice -n 40 make all JOBSNUM=16 DECI=6   MPFR=0 SSE=0 FETS=-8 THE_FAST=EXPERIMENTALLY_FAST=OFF    && echo "\e[96m==intel=A 1 ===\e[0m"  &&  ./examples/yade --test ; ./examples/yade --check  ; ./examples/yade --quickperformance -j16  ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
make uninstall && nice -n 40 make all JOBSNUM=16 DECI=6   MPFR=0 SSE=0 FETS=-8 THE_FAST=DANGEROUSLY_FAST=OFF SSE=1 && echo "\e[96m==intel=A 2 ===\e[0m"  &&  ./examples/yade --test ; ./examples/yade --check  ; ./examples/yade --quickperformance -j16  ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
make uninstall && nice -n 40 make all JOBSNUM=16 DECI=6   MPFR=0 SSE=0 FETS=-8 THE_FAST=EXPERIMENTALLY_FAST=ON     && echo "\e[96m==intel=A 3 ===\e[0m"  &&  ./examples/yade --test ; ./examples/yade --check  ; ./examples/yade --quickperformance -j16  ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
make uninstall && nice -n 40 make all JOBSNUM=16 DECI=6   MPFR=0 SSE=0 FETS=-8 THE_FAST=DANGEROUSLY_FAST=ON        && echo "\e[96m==intel=A 4 ===\e[0m"  &&  ./examples/yade --test ; ./examples/yade --check  ; ./examples/yade --quickperformance -j16  ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
make uninstall && nice -n 40 make all JOBSNUM=16 DECI=6   MPFR=0 SSE=0 FETS=-8 THE_FAST=DANGEROUSLY_FAST=ON  SSE=1 && echo "\e[96m==intel=A 5 ===\e[0m"  &&  ./examples/yade --test ; ./examples/yade --check  ; ./examples/yade --quickperformance -j16  ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
date
make uninstall && nice -n 40 make all JOBSNUM=16 DECI=15  MPFR=0 SSE=0 FETS=-8 THE_FAST=EXPERIMENTALLY_FAST=OFF    && echo "\e[96m==intel=B 1 ===\e[0m"  &&  ./examples/yade --test ; ./examples/yade --check  ; ./examples/yade --quickperformance -j16  ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
make uninstall && nice -n 40 make all JOBSNUM=16 DECI=15  MPFR=0 SSE=0 FETS=-8 THE_FAST=DANGEROUSLY_FAST=OFF SSE=1 && echo "\e[96m==intel=B 2 ===\e[0m"  &&  ./examples/yade --test ; ./examples/yade --check  ; ./examples/yade --quickperformance -j16  ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
make uninstall && nice -n 40 make all JOBSNUM=16 DECI=15  MPFR=0 SSE=0 FETS=-8 THE_FAST=EXPERIMENTALLY_FAST=ON     && echo "\e[96m==intel=B 3 ===\e[0m"  &&  ./examples/yade --test ; ./examples/yade --check  ; ./examples/yade --quickperformance -j16  ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
make uninstall && nice -n 40 make all JOBSNUM=16 DECI=15  MPFR=0 SSE=0 FETS=-8 THE_FAST=DANGEROUSLY_FAST=ON        && echo "\e[96m==intel=B 4 ===\e[0m"  &&  ./examples/yade --test ; ./examples/yade --check  ; ./examples/yade --quickperformance -j16  ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
make uninstall && nice -n 40 make all JOBSNUM=16 DECI=15  MPFR=0 SSE=0 FETS=-8 THE_FAST=DANGEROUSLY_FAST=ON  SSE=1 && echo "\e[96m==intel=B 5 ===\e[0m"  &&  ./examples/yade --test ; ./examples/yade --check  ; ./examples/yade --quickperformance -j16  ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
date
make uninstall && nice -n 40 make all JOBSNUM=16 DECI=18  MPFR=0 SSE=0 FETS=-8 THE_FAST=EXPERIMENTALLY_FAST=OFF    && echo "\e[96m==intel=C 1 ===\e[0m"  &&  ./examples/yade --test ; ./examples/yade --check  ; ./examples/yade --quickperformance -j16  ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
make uninstall && nice -n 40 make all JOBSNUM=16 DECI=18  MPFR=0 SSE=0 FETS=-8 THE_FAST=DANGEROUSLY_FAST=OFF SSE=1 && echo "\e[96m==intel=C 2 ===\e[0m"  &&  ./examples/yade --test ; ./examples/yade --check  ; ./examples/yade --quickperformance -j16  ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
make uninstall && nice -n 40 make all JOBSNUM=16 DECI=18  MPFR=0 SSE=0 FETS=-8 THE_FAST=EXPERIMENTALLY_FAST=ON     && echo "\e[96m==intel=C 3 ===\e[0m"  &&  ./examples/yade --test ; ./examples/yade --check  ; ./examples/yade --quickperformance -j16  ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
make uninstall && nice -n 40 make all JOBSNUM=16 DECI=18  MPFR=0 SSE=0 FETS=-8 THE_FAST=DANGEROUSLY_FAST=ON        && echo "\e[96m==intel=C 4 ===\e[0m"  &&  ./examples/yade --test ; ./examples/yade --check  ; ./examples/yade --quickperformance -j16  ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
make uninstall && nice -n 40 make all JOBSNUM=16 DECI=18  MPFR=0 SSE=0 FETS=-8 THE_FAST=DANGEROUSLY_FAST=ON  SSE=1 && echo "\e[96m==intel=C 5 ===\e[0m"  &&  ./examples/yade --test ; ./examples/yade --check  ; ./examples/yade --quickperformance -j16  ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
date
make uninstall && nice -n 40 make all JOBSNUM=16 DECI=33  MPFR=0 SSE=0 FETS=-8 THE_FAST=EXPERIMENTALLY_FAST=OFF    && echo "\e[96m==intel=D 1 ===\e[0m"  &&  ./examples/yade --test ; ./examples/yade --check  ; ./examples/yade --quickperformance -j16  ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
make uninstall && nice -n 40 make all JOBSNUM=16 DECI=33  MPFR=0 SSE=0 FETS=-8 THE_FAST=DANGEROUSLY_FAST=OFF SSE=1 && echo "\e[96m==intel=D 2 ===\e[0m"  &&  ./examples/yade --test ; ./examples/yade --check  ; ./examples/yade --quickperformance -j16  ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
make uninstall && nice -n 40 make all JOBSNUM=16 DECI=33  MPFR=0 SSE=0 FETS=-8 THE_FAST=EXPERIMENTALLY_FAST=ON     && echo "\e[96m==intel=D 3 ===\e[0m"  &&  ./examples/yade --test ; ./examples/yade --check  ; ./examples/yade --quickperformance -j16  ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
make uninstall && nice -n 40 make all JOBSNUM=16 DECI=33  MPFR=0 SSE=0 FETS=-8 THE_FAST=DANGEROUSLY_FAST=ON        && echo "\e[96m==intel=D 4 ===\e[0m"  &&  ./examples/yade --test ; ./examples/yade --check  ; ./examples/yade --quickperformance -j16  ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
make uninstall && nice -n 40 make all JOBSNUM=16 DECI=33  MPFR=0 SSE=0 FETS=-8 THE_FAST=DANGEROUSLY_FAST=ON  SSE=1 && echo "\e[96m==intel=D 5 ===\e[0m"  &&  ./examples/yade --test ; ./examples/yade --check  ; ./examples/yade --quickperformance -j16  ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
date
make uninstall && nice -n 40 make all JOBSNUM=16 DECI=40  MPFR=1 SSE=0 FETS=-8 THE_FAST=EXPERIMENTALLY_FAST=OFF    && echo "\e[96m==intel=E 1 ===\e[0m"  &&  ./examples/yade --test ; ./examples/yade --check  ; ./examples/yade --quickperformance -j16  ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
make uninstall && nice -n 40 make all JOBSNUM=16 DECI=40  MPFR=1 SSE=0 FETS=-8 THE_FAST=DANGEROUSLY_FAST=OFF SSE=1 && echo "\e[96m==intel=E 2 ===\e[0m"  &&  ./examples/yade --test ; ./examples/yade --check  ; ./examples/yade --quickperformance -j16  ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
make uninstall && nice -n 40 make all JOBSNUM=16 DECI=40  MPFR=1 SSE=0 FETS=-8 THE_FAST=EXPERIMENTALLY_FAST=ON     && echo "\e[96m==intel=E 3 ===\e[0m"  &&  ./examples/yade --test ; ./examples/yade --check  ; ./examples/yade --quickperformance -j16  ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
make uninstall && nice -n 40 make all JOBSNUM=16 DECI=40  MPFR=1 SSE=0 FETS=-8 THE_FAST=DANGEROUSLY_FAST=ON        && echo "\e[96m==intel=E 4 ===\e[0m"  &&  ./examples/yade --test ; ./examples/yade --check  ; ./examples/yade --quickperformance -j16  ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
make uninstall && nice -n 40 make all JOBSNUM=16 DECI=40  MPFR=1 SSE=0 FETS=-8 THE_FAST=DANGEROUSLY_FAST=ON  SSE=1 && echo "\e[96m==intel=E 5 ===\e[0m"  &&  ./examples/yade --test ; ./examples/yade --check  ; ./examples/yade --quickperformance -j16  ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
date
make uninstall && nice -n 40 make all JOBSNUM=16 DECI=150 MPFR=1 SSE=0 FETS=-8 THE_FAST=EXPERIMENTALLY_FAST=OFF    && echo "\e[96m==intel=F 1 ===\e[0m"  &&  ./examples/yade --test ; ./examples/yade --check  ; ./examples/yade --quickperformance -j16  ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
make uninstall && nice -n 40 make all JOBSNUM=16 DECI=150 MPFR=1 SSE=0 FETS=-8 THE_FAST=DANGEROUSLY_FAST=OFF SSE=1 && echo "\e[96m==intel=F 2 ===\e[0m"  &&  ./examples/yade --test ; ./examples/yade --check  ; ./examples/yade --quickperformance -j16  ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
make uninstall && nice -n 40 make all JOBSNUM=16 DECI=150 MPFR=1 SSE=0 FETS=-8 THE_FAST=EXPERIMENTALLY_FAST=ON     && echo "\e[96m==intel=F 3 ===\e[0m"  &&  ./examples/yade --test ; ./examples/yade --check  ; ./examples/yade --quickperformance -j16  ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
make uninstall && nice -n 40 make all JOBSNUM=16 DECI=150 MPFR=1 SSE=0 FETS=-8 THE_FAST=DANGEROUSLY_FAST=ON        && echo "\e[96m==intel=F 4 ===\e[0m"  &&  ./examples/yade --test ; ./examples/yade --check  ; ./examples/yade --quickperformance -j16  ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
make uninstall && nice -n 40 make all JOBSNUM=16 DECI=150 MPFR=1 SSE=0 FETS=-8 THE_FAST=DANGEROUSLY_FAST=ON  SSE=1 && echo "\e[96m==intel=F 5 ===\e[0m"  &&  ./examples/yade --test ; ./examples/yade --check  ; ./examples/yade --quickperformance -j16  ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
date
make uninstall && nice -n 40 make all JOBSNUM=16 DECI=40  MPFR=0 SSE=0 FETS=-8 THE_FAST=EXPERIMENTALLY_FAST=OFF    && echo "\e[96m==intel=G 1 ===\e[0m"  &&  ./examples/yade --test ; ./examples/yade --check  ; ./examples/yade --quickperformance -j16  ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
make uninstall && nice -n 40 make all JOBSNUM=16 DECI=40  MPFR=0 SSE=0 FETS=-8 THE_FAST=DANGEROUSLY_FAST=OFF SSE=1 && echo "\e[96m==intel=G 2 ===\e[0m"  &&  ./examples/yade --test ; ./examples/yade --check  ; ./examples/yade --quickperformance -j16  ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
make uninstall && nice -n 40 make all JOBSNUM=16 DECI=40  MPFR=0 SSE=0 FETS=-8 THE_FAST=EXPERIMENTALLY_FAST=ON     && echo "\e[96m==intel=G 3 ===\e[0m"  &&  ./examples/yade --test ; ./examples/yade --check  ; ./examples/yade --quickperformance -j16  ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
make uninstall && nice -n 40 make all JOBSNUM=16 DECI=40  MPFR=0 SSE=0 FETS=-8 THE_FAST=DANGEROUSLY_FAST=ON        && echo "\e[96m==intel=G 4 ===\e[0m"  &&  ./examples/yade --test ; ./examples/yade --check  ; ./examples/yade --quickperformance -j16  ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
make uninstall && nice -n 40 make all JOBSNUM=16 DECI=40  MPFR=0 SSE=0 FETS=-8 THE_FAST=DANGEROUSLY_FAST=ON  SSE=1 && echo "\e[96m==intel=G 5 ===\e[0m"  &&  ./examples/yade --test ; ./examples/yade --check  ; ./examples/yade --quickperformance -j16  ; ./examples/yade --stdperformance -j16 ; echo "\n\n\n\n=====\n\n\n" 
date

