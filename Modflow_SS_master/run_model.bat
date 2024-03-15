
@echo off > nul

del Gondebeek_SS.hds 

del Gondebeek_SS.cbc 

mf6.exe Gondebeek_SS.nam 

mf6mod2obs < simulated.in 
