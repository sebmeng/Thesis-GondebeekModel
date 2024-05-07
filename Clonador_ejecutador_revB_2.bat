::Autores: GDM (JCortinez@Arcadis2016)
::	   GDM (GLetelier@Arcadis2016)
:: revB (versión con RoBOCOPY)
:: Clonador ejecutador
:: Basicamente el script clona una carpeta_base y la replica en el directorio donde
:: se encuentra el script, dado un numero inicial y final de las realizaciones
:: Posteriormente corre un archivo cualquiera en un nuevo terminal, en este caso 
:: un script de un modelo Modflow

@echo off 
color 17

::set archivo_ejecutable=4655_RunModel.bat
::set archivo_ejecutable=USGs_1 4655_Sim_02
set archivo_ejecutable=RunSlave.bat

set carpeta_base=0
echo Ruta de la carpeta original:
set parent_folder=%~dp0%carpeta_base%
echo %parent_folder%

echo Ruta de salida:
set output_folder=%~dp0
echo %output_folder%

echo Numero INICIAL de carpetas replicadas:
set N_Folderi=
set /P N_Folderi=Type input: %=%

echo Numero FINAL de carpetas replicadas:
set N_Folderf=
set /P N_Folderf=Type input: %=%

if exist "%parent_folder%" goto validador
ECHO ERROR: Carpeta no encontrada

:validador
if exist "%output_folder%" goto copia
ECHO ERROR: Carpeta de salida no encotrada

:copia
ECHO Carpeta principal encontrada...iniciando la copia
setlocal EnableDelayedExpansion
FOR /L %%A IN (%N_Folderi%,1,%N_Folderf%) DO (
	robocopy %parent_folder% %output_folder%%\%%%A% /E
	cd %output_folder%%\%%%A%
	start %output_folder%%\%%%A% %\%archivo_ejecutable%
	echo %output_folder%%\%%%A% %\%archivo_ejecutable%
)
echo ===================================
ECHO ======  Fin del script   ==========
echo ===================================
pause > nul