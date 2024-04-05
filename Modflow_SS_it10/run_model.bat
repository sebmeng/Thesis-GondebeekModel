
@echo off > nul

del Gondebeek_SS.hds > nul

del Gondebeek_SS.cbc > nul

mf6.exe Gondebeek_SS.nam > nul

mf6mod2obs < simulated.in > nul
