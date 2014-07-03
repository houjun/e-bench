#include <stdio.h>
#include <stdlib.h>
#include "hdf5.h"

#define DATAFILE    "/scratch2/scratchdirs/houhun/data/brain/EC6_CV.h5"
#define ROOTPROC    0
#define MAXDIM      3
#define LABEL_LEN   6

typedef struct ECoGMeta {

    // Dimension
    int     ECoGIndx_ndims;
    int     EIndx_ndims;   

    // Data size of each dimension
    hsize_t ECoGIndx_dim[MAXDIM];       // Offset range of each trial
    hsize_t EIndx_dim[MAXDIM];          // Type ID of each trial
    hsize_t ELbls_dim[MAXDIM];          // Label name of each type ID

    // actual metadata
    double* ECoGIndx_data;
    double* EIndx_data;
    char*  ELbls_data;

} ECoGMeta;

herr_t root_get_metadata(hid_t file_id, ECoGMeta* metadata) 
{

    int     i;
    hid_t   ECoGIndx_id, EIndx_id, ELbls_id;
    hid_t   ECoGIndx_space, EIndx_space, ELbls_space;
    hsize_t ECoGIndx_size=1, EIndx_size=1, ELbls_size=1;
    herr_t      status;
     

    // Open datasets, assuming the file is already open and file_id is correct
    ECoGIndx_id = H5Dopen(file_id, "/Descriptors/Event_ECoGIndx", H5P_DEFAULT);
    EIndx_id    = H5Dopen(file_id, "/Descriptors/Event_EIndx", H5P_DEFAULT);
    ELbls_id    = H5Dopen(file_id, "/Descriptors/Event_ELbls", H5P_DEFAULT);

    // Get space
    ECoGIndx_space = H5Dget_space(ECoGIndx_id);
    EIndx_space    = H5Dget_space(EIndx_id   );
    ELbls_space    = H5Dget_space(ELbls_id   );

    // Get dims of each datasets
    metadata->ECoGIndx_ndims = H5Sget_simple_extent_ndims(ECoGIndx_space);
    metadata->EIndx_ndims    = H5Sget_simple_extent_ndims(EIndx_space);

    // Get data size
    H5Sget_simple_extent_dims(ECoGIndx_space, metadata->ECoGIndx_dim, NULL);
    H5Sget_simple_extent_dims(EIndx_space,    metadata->EIndx_dim, NULL);
    H5Sget_simple_extent_dims(ELbls_space,    metadata->ELbls_dim, NULL);

    // Calculate how much space we need
    for (i = 0; i < metadata->ECoGIndx_ndims; i++) 
        ECoGIndx_size *= metadata->ECoGIndx_dim[i];

    for (i = 0; i < metadata->EIndx_ndims; i++) 
        EIndx_size    *= metadata->ECoGIndx_dim[i];

    // Allocate memory
    metadata->ECoGIndx_data = (double*)malloc(ECoGIndx_size*sizeof(double));
    metadata->EIndx_data    = (double*)malloc(EIndx_size*sizeof(double));
    metadata->ELbls_data    = (char*)malloc(metadata->ELbls_dim[0]*LABEL_LEN*sizeof(double));

    // Read all metadata
    H5Dread(ECoGIndx_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, metadata->ECoGIndx_data);
    H5Dread(EIndx_id   , H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, metadata->EIndx_data);
//    H5Dread(ELbls_id   , H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, metadata->ELbls_data);


    // Close
    status = H5Dclose(ECoGIndx_id);
    status = H5Dclose(EIndx_id   );
    status = H5Dclose(ELbls_id   );

    return status;

}

int main(int argc, char* argv[]) {

    int         proc_num, my_rank;
    int         i, j, k;
    
    ECoGMeta*   metadata;
    hid_t       ECoGData_id, file_id;  
    herr_t      status;

    // Start MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    // All processes open the data file.
    file_id = H5Fopen(DATAFILE, H5F_ACC_RDWR, H5P_DEFAULT);

    // Only root process open and read all index datasets
    if (my_rank == ROOTPROC) { 
     
        metadata = malloc(sizeof(ECoGMeta));
        status = root_get_metadata(file_id, metadata);

        printf("Event_ECoGIndx: [0]:%.0f [1]:%.0f\n", metadata->ECoGIndx_data[0], metadata->ECoGIndx_data[1]);
        printf("Event_ECoGIndx: [3510]:%.0f [3511]:%.0f\n", metadata->ECoGIndx_data[3510], metadata->ECoGIndx_data[3511]);
        printf("EIndx_data: [0]:%.0f [1755]:%.0f\n", metadata->EIndx_data[0], metadata->EIndx_data[1755]);

    }



 
    // Open an existing dataset.
    ECoGData_id = H5Dopen(file_id, "/Data/ECoG", H5P_DEFAULT);
 
 
    // Close everything. 
    status = H5Dclose(ECoGData_id);
    status = H5Fclose(file_id);


    H5close();
    MPI_Finalize();
    return 0;

}

