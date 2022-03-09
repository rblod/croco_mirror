#-------------------------------------------------------------------------------
# OCE
#-------------------------------------------------------------------------------
# namelist [Info: grid size is directly read in oce_compile.sh and cpl_nam.sh ]

# Online Compilation
# Creates croco executable depending on this namelist. 
#     In cppdef.h options that can be edit are OA_COUPLING/OW_COUPLING/BULK_FLUX/ONLINE/CLIMATOLOGY/FRC_BRY/TIDES/XIOS/AGRIF/AGRIF_2WAY
#     In param.h modify the size grid, edit the number of proc in x and y direction with intel given in myjob.sh
export ONLINE_COMP=0

# Time steps
export DT_OCE=600
export NDTFAST=60

# Parameter 
export hmin=75; # minimum water depth in CROCO, delimiting coastline in WW3 

# domains
export AGRIFZ=0
export AGRIF_2WAY="FALSE"

# MPI_NOLAND
export MPI_NOLAND="FALSE"
export MY_NODES=18 # ONLY if MPI_NOLAND is "TRUE". It replaces NP_OCEX*NP_OCEY

# forcing files
export ini_ext='ini_SODA' # ini extension file (ini_SODA,...)
export bdy_ext='bry_SODA' # bry extension file (clm_SODA,bry_SODA,...)
export surfrc_flag="TRUE" # Flag if surface forcing is needed (FALSE if coupling with the atmosphere, TRUE otherwise)
export interponline=0 # switch (1=on, 0=off) for online surface interpolation. Only works with MONTHLY input files!
export frc_ext='blk_CFSR' # surface forcing extension(blk_CFSR, frc_CFSR,...). If interponline=1 precise the type (ECMWF, CFSR,...). If AROME/ARPEGE, specifies the file name AROME_*.nc (Inthis case, by default we define 24 records per days)
export tide_flag="FALSE" # the forcing extension must be blk_??? otherwise tide forcing overwrites it 
export river_flag="FALSE"

# output settings
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#                                          WARNING                                       ! 
# When XIOS is activated the following values (for the model) are not taken into account !
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
export oce_his_sec=3600     # history output interval (in number of second) 
export oce_avg_sec=3600     # average output interval (in number of second) 


