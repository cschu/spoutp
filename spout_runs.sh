#!/bin/bash
QUEUE=TSL-Test256 #128
SPOUTP=/usr/users/sl/schudomc/spoutp/spoutp_cluster #_debug


#bsub -q $QUEUE "$SPOUTP /tsl/services/galaxy/dist/galaxy-dist/database/files/023/dataset_23932.dat Go5_spot_24h_C  10"
#bsub -q $QUEUE "$SPOUTP /tsl/services/galaxy/dist/galaxy-dist/database/files/023/dataset_23823.dat classified_48h_G05_Spray 10"
#bsub -q $QUEUE "$SPOUTP /tsl/services/galaxy/dist/galaxy-dist/database/files/023/dataset_23826.dat unclassified_48h_G05_Spray 10"
#exit

bsub -q $QUEUE "$SPOUTP /tsl/services/galaxy/dist/galaxy-dist/database/files/023/dataset_23759.dat classified_0h_G05_Spray 10; $SPOUTP /tsl/services/galaxy/dist/galaxy-dist/database/files/023/dataset_23761.dat unclassified_0h_G05_Spray 10; $SPOUTP /tsl/services/galaxy/dist/galaxy-dist/database/files/023/dataset_23763.dat classified_24h_G05_Spray 10; $SPOUTP /tsl/services/galaxy/dist/galaxy-dist/database/files/023/dataset_23765.dat unclassified_24h_G05_Spray 10; $SPOUTP /tsl/services/galaxy/dist/galaxy-dist/database/files/023/dataset_23776.dat classified_36h_G05_Spray 10; $SPOUTP /tsl/services/galaxy/dist/galaxy-dist/database/files/023/dataset_23778.dat unclassified_36h_G05_Spray 10; $SPOUTP /tsl/services/galaxy/dist/galaxy-dist/database/files/023/dataset_23823.dat classified_48h_G05_Spray 10;  $SPOUTP /tsl/services/galaxy/dist/galaxy-dist/database/files/023/dataset_23826.dat unclassified_48h_G05_Spray 10"

bsub -q $QUEUE "$SPOUTP /tsl/services/galaxy/dist/galaxy-dist/database/files/023/dataset_23828.dat classified_72h_G05_Spray 10; $SPOUTP /tsl/services/galaxy/dist/galaxy-dist/database/files/023/dataset_23830.dat unclassified_72h_G05_Spray 10; $SPOUTP /tsl/services/galaxy/dist/galaxy-dist/database/files/023/dataset_23846.dat classified_12h_G05_Spray 10; $SPOUTP /tsl/services/galaxy/dist/galaxy-dist/database/files/023/dataset_23848.dat unclassified_12h_G05_Spray 10; $SPOUTP /tsl/services/galaxy/dist/galaxy-dist/database/files/023/dataset_23687.dat classified_96+168h_G05_Spray 10; $SPOUTP /tsl/services/galaxy/dist/galaxy-dist/database/files/023/dataset_23689.dat unclassified_96+168h_G05_Spray 10; $SPOUTP /tsl/services/galaxy/dist/galaxy-dist/database/files/022/dataset_22378.dat classified_mi168 10; $SPOUTP /tsl/services/galaxy/dist/galaxy-dist/database/files/022/dataset_22380.dat unclassified_mi168 10; $SPOUTP /tsl/services/galaxy/dist/galaxy-dist/database/files/022/dataset_22359.dat unclassified_mi12 10"

#bsub -q $QUEUE "$SPOUTP /tsl/services/galaxy/dist/galaxy-dist/database/files/022/dataset_22359.dat unclassified_mi12 10" 
#bsub -q $QUEUE "$SPOUTP /tsl/services/galaxy/dist/galaxy-dist/database/files/022/dataset_22357.dat classified_mi12 10" 

 bsub -q $QUEUE "$SPOUTP /tsl/services/galaxy/dist/galaxy-dist/database/files/022/dataset_22373.dat classified_mi72 10; $SPOUTP /tsl/services/galaxy/dist/galaxy-dist/database/files/022/dataset_22375.dat unclassified_mi72 10; $SPOUTP /tsl/services/galaxy/dist/galaxy-dist/database/files/022/dataset_22369.dat classified_mi36 10; $SPOUTP /tsl/services/galaxy/dist/galaxy-dist/database/files/022/dataset_22371.dat unclassified_mi36 10; $SPOUTP /tsl/services/galaxy/dist/galaxy-dist/database/files/022/dataset_22365.dat classified_mi24 10; $SPOUTP /tsl/services/galaxy/dist/galaxy-dist/database/files/022/dataset_22367.dat unclassified_mi24 10; $SPOUTP /tsl/services/galaxy/dist/galaxy-dist/database/files/022/dataset_22351.dat classified_mi0 10; $SPOUTP /tsl/services/galaxy/dist/galaxy-dist/database/files/022/dataset_22353.dat unclassified_mi0 10; $SPOUTP /tsl/services/galaxy/dist/galaxy-dist/database/files/022/dataset_22357.dat classified_mi12 10"

bsub -q $QUEUE "$SPOUTP /tsl/services/galaxy/dist/galaxy-dist/database/files/023/dataset_23928.dat Go5_spot_0h_C   10; $SPOUTP /tsl/services/galaxy/dist/galaxy-dist/database/files/023/dataset_23930.dat Go5_spot_0h_U   10; $SPOUTP /tsl/services/galaxy/dist/galaxy-dist/database/files/023/dataset_23932.dat Go5_spot_24h_C  10; $SPOUTP /tsl/services/galaxy/dist/galaxy-dist/database/files/023/dataset_23934.dat Go5_spot_24h_U  10; $SPOUTP /tsl/services/galaxy/dist/galaxy-dist/database/files/023/dataset_23936.dat Go5_spot_48h_C  10; $SPOUTP /tsl/services/galaxy/dist/galaxy-dist/database/files/023/dataset_23938.dat Go5_spot_48h_U  10; $SPOUTP /tsl/services/galaxy/dist/galaxy-dist/database/files/023/dataset_23940.dat Go5_spot_468_C 10; $SPOUTP /tsl/services/galaxy/dist/galaxy-dist/database/files/023/dataset_23942.dat Go5_spot_468_U 10"


#bsub -q $QUEUE "$SPOUTP /tsl/services/galaxy/dist/galaxy-dist/database/files/022/dataset_22378.dat classified_mi168 10"
#bsub -q $QUEUE "$SPOUTP /tsl/services/galaxy/dist/galaxy-dist/database/files/022/dataset_22380.dat unclassified_mi168 10"







