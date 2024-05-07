
@echo off > nul

del Gondebeek.hds > nul

del Gondebeek.cbc > nul

mf6.exe Gondebeek.nam 

mf6mod2obs < simulated_PEST.in > nul
