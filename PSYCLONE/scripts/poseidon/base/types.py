##########################################################
#  CROCO PSYCLONE scripts, under CeCILL-C
#  From Sébastien Valat (INRIA) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
# python
from enum import Enum

##########################################################
class AccessMode(Enum):
    UNDEFINED = 0
    READ = 1
    WRITE = 2
    READ_WRITE = 3
