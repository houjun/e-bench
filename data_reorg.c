#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
#include "hdf5.h"
#include "e-bench.h"

int read_metadata(char* filename, ECoGMeta* metadata)
{
    int     i;
    hid_t   file_id;
    hid_t   ECoGIndx_id, EIndx_id, ELbls_id, ELbls_memtype;
    hid_t   ECoGIndx_space, EIndx_space, ELbls_space;
    herr_t      status;

    file_id = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);

    // Open datasets, assuming the file is already open and file_id is correct
    ECoGIndx_id = H5Dopen(file_id, "/Descriptors/Event_ECoGIndx", H5P_DEFAULT);
    EIndx_id    = H5Dopen(file_id, "/Descriptors/Event_EIndx", H5P_DEFAULT);
    ELbls_id    = H5Dopen(file_id, "/Descriptors/Event_ELbls", H5P_DEFAULT);

    // Get space
    ECoGIndx_space = H5Dget_space(ECoGIndx_id);
    EIndx_space    = H5Dget_space(EIndx_id   );
    ELbls_space    = H5Dget_space(ELbls_id   );

    // Get data size and dim info
    metadata->ECoGIndx_ndims = H5Sget_simple_extent_dims(ECoGIndx_space, metadata->ECoGIndx_dim, NULL);
    metadata->EIndx_ndims    = H5Sget_simple_extent_dims(EIndx_space,    metadata->EIndx_dim, NULL);
    metadata->ELbls_ndims    = H5Sget_simple_extent_dims(ELbls_space,    metadata->ELbls_dim, NULL);

    metadata->ECoGIndx_size  = 1;
    // Calculate how much space we need
    for (i = 0; i < metadata->ECoGIndx_ndims; i++)
        metadata->ECoGIndx_size *= metadata->ECoGIndx_dim[i];

    metadata->EIndx_size    = metadata->EIndx_dim[0];

    // Allocate memory
    metadata->ECoGIndx_data = (double*)malloc(metadata->ECoGIndx_size*sizeof(double));
    metadata->EIndx_data    = (double*)malloc(metadata->EIndx_size*sizeof(double));

    char** tmpELbls_data    = (char**)malloc(metadata->ELbls_dim[0]*sizeof(char*));

    metadata->ELbls_data    = (char**)malloc(metadata->ELbls_dim[0]*sizeof(char*));
    for (i = 0; i < metadata->ELbls_dim[0]; i++) {
        metadata->ELbls_data[i] = (char*)malloc(LABEL_LEN*sizeof(char));
    }

    // Read all metadata
    status = H5Dread(ECoGIndx_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, metadata->ECoGIndx_data);
    status = H5Dread(EIndx_id   , H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, metadata->EIndx_data);

    // Special case for reading label strings.
    ELbls_memtype = H5Tcopy(H5T_C_S1);
    status = H5Tset_size(ELbls_memtype, H5T_VARIABLE);
    status = H5Dread(ELbls_id, ELbls_memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmpELbls_data);
    for (i = 0; i < metadata->ELbls_dim[0]; i++) {
        strcpy(metadata->ELbls_data[i], tmpELbls_data[i]);
    }


    // Close
    status = H5Dvlen_reclaim (ELbls_memtype, ELbls_space, H5P_DEFAULT, tmpELbls_data);
    free(tmpELbls_data);
    status = H5Dclose(ECoGIndx_id);
    status = H5Dclose(EIndx_id);
    status = H5Dclose(ELbls_id);
    status = H5Sclose(ECoGIndx_space);
    status = H5Sclose(EIndx_space);
    status = H5Sclose(ELbls_space);
    status = H5Tclose(ELbls_memtype);
    status = H5Fclose(file_id);
}

int data_reorg(char* filename, ECoGMeta* metadata)
{
    herr_t      status;
    hid_t       file_id, ECoGData_id, ECoGData_space, ECoGData_memspace, file_space;
    char        new_filename[HNAME_MAX];
    hid_t       new_file_id, new_ECoGData_id, new_ECoGData_g, new_ECoGData_space, new_ECoGData_memspace;
    hsize_t     i, j, k;
    double*     new_EIndx = (double*)malloc(metadata->EIndx_size*sizeof(double));
    double*     new_ECoGIndx = (double*)malloc(metadata->ECoGIndx_size*sizeof(double));

    hsize_t     ECoGData_dim[MAXDIM], file_sel[MAXDIM];
    hsize_t     my_count[MAXDIM], my_offset[MAXDIM];

    sprintf(new_filename,"%s_new", filename);
    printf("%s\n", new_filename);
    // Open file, dataset
    file_id              = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    ECoGData_id          = H5Dopen(file_id, DSET_NAME, H5P_DEFAULT);
    ECoGData_space       = H5Dget_space(ECoGData_id);


    H5Sget_simple_extent_dims(ECoGData_space, ECoGData_dim, NULL);

    my_offset[1]         = 0; 
    my_count[0]          = metadata->ECoGIndx_data[1] - metadata->ECoGIndx_data[0] + 1;
    my_count[1]          = ECoGData_dim[1]; 

    // Open new file, dataset
    new_file_id          = H5Fcreate(new_filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    new_ECoGData_g       = H5Gcreate(new_file_id, "/Data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    new_ECoGData_space   = H5Screate_simple(2, ECoGData_dim, NULL);
    new_ECoGData_id      = H5Dcreate(new_ECoGData_g, "/Data/ECoG", H5T_IEEE_F64LE, new_ECoGData_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);


    // Allocate space for each trial
    hsize_t trial_size   = my_count[0] * my_count[1];
    double* ECoGData     = (double*)malloc(sizeof(double)*trial_size);

    printf("trial size:%llu EIndx_dim[0]=%llu\n", trial_size, metadata->EIndx_dim[0]);

    ECoGData_memspace    = H5Screate_simple(2, my_count, NULL);

    hsize_t write_offset = 0;
    // for all trial type
    k = 0;
    hsize_t tmp;
    for (j = 0; j < metadata->ELbls_dim[0]; j++) {
   
        // for all trials
        for (i = 0; i < metadata->EIndx_dim[0]; i++) {
            tmp = (hsize_t)(metadata->EIndx_data[i]) - 1;
            if ( tmp == j ) {

                my_count[0]  = metadata->ECoGIndx_data[i*2+1] - metadata->ECoGIndx_data[i*2] + 1;
                my_offset[0] = (metadata->ECoGIndx_data[i*2] - 1);

                if (my_count[0] != 301) {
                    printf("%d: %f - %f = %llu\n", i, (metadata->ECoGIndx_data[i*2+1]), (metadata->ECoGIndx_data[i*2]), my_count[0]);
                }
                //printf("Read - Offsets: [0]:%llu [1]:%llu\n", my_offset[0],my_offset[1]);
                //printf("Read - Counts:  [0]:%llu [1]:%llu\n", my_count[0],my_count[1]);

                // read the data
                status       = H5Sselect_hyperslab(ECoGData_space, H5S_SELECT_SET, my_offset, NULL, my_count, NULL);
                status       = H5Dread(ECoGData_id, H5T_IEEE_F64LE, ECoGData_memspace, ECoGData_space, H5P_DEFAULT, ECoGData);
                //printf("[0]:%f,[1]:%f\n",ECoGData[0],ECoGData[1]);

                // write to new place
                my_offset[0] = write_offset;
                file_space   = H5Dget_space(new_ECoGData_id);
                status       = H5Sselect_hyperslab(file_space, H5S_SELECT_SET, my_offset, NULL, my_count, NULL);
                status       = H5Dwrite(new_ECoGData_id, H5T_IEEE_F64LE, ECoGData_memspace, file_space, H5P_DEFAULT, ECoGData);

                //printf("Write - Offsets: [0]:%llu [1]:%llu\n", my_offset[0],my_offset[1]);
                //printf("Write - Counts:  [0]:%llu [1]:%llu\n", my_count[0],my_count[1]);

                // update index
                new_ECoGIndx[k*2]   = write_offset + 1;
                new_ECoGIndx[k*2+1] = write_offset + my_count[0];
                new_EIndx[k++]      = j+1;

                write_offset += my_count[0];
            }
            
        }

    }

    printf("k=%d\n",k);

    new_ECoGData_g       = H5Gcreate(new_file_id, "/Descriptors", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    new_ECoGData_space   = H5Screate_simple(2, metadata->ECoGIndx_dim, NULL);
    new_ECoGData_id      = H5Dcreate(new_ECoGData_g, "/Descriptors/Event_ECoGIndx", H5T_IEEE_F64LE, new_ECoGData_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status               = H5Dwrite(new_ECoGData_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, new_ECoGIndx);

    new_ECoGData_space   = H5Screate_simple(2, metadata->EIndx_dim, NULL);
    new_ECoGData_id      = H5Dcreate(new_ECoGData_g, "/Descriptors/Event_EIndx", H5T_IEEE_F64LE, new_ECoGData_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status               = H5Dwrite(new_ECoGData_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, new_EIndx);

    hid_t       filetype = H5Tcopy (H5T_C_S1);
    status               = H5Tset_size (filetype, H5T_VARIABLE);
    hid_t memtype        = H5Tcopy (H5T_C_S1);
    status               = H5Tset_size (memtype, H5T_VARIABLE);

    hsize_t dims[1];
    dims[0] = metadata->ELbls_dim[0];
    hid_t space          = H5Screate_simple (1, dims, NULL);
    hid_t dset           = H5Dcreate (new_file_id, "/Descriptors/Event_ELbls", filetype, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status               = H5Dwrite (dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, metadata->ELbls_data);


    status = H5Dclose(new_ECoGData_id);
    status = H5Sclose(new_ECoGData_space);
    status = H5Fclose(new_file_id);

    status = H5Dclose(ECoGData_id);
    status = H5Sclose(ECoGData_space);
    status = H5Sclose(ECoGData_memspace);
    status = H5Fclose(file_id);


}

int main(int argc, char* argv[])
{
    ECoGMeta* metadata = (ECoGMeta*)malloc(sizeof(ECoGMeta));

    char* fname;

    if(argc == 1)
        fname = "EC6_CV.h5";
    else
        fname = argv[1];
    read_metadata(fname, metadata);

    data_reorg(fname, metadata);

    return 0;
}

