#!/bin/csh

# file to set up csh/tcsh shell environment to use R

source /usr/local/apps/java/java18.csh

if !($?PATH) then
    setenv PATH /usr/local/apps/R/em64t/R-3.4.2/bin
else
    setenv PATH /usr/local/apps/R/em64t/R-3.4.2/bin:$PATH
endif

