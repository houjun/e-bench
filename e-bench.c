#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
#include "hdf5.h"
#include "e-bench.h"

int main(int argc, char* argv[]) {

    
    herr_t      status;
    hid_t       file_id, ECoGData_id, ECoGData_space, ECoGData_memspace, plist_id;
    hsize_t     i, j, ECoGData_dim[MAXDIM], my_count[MAXDIM], my_offset[MAXDIM], my_trials, total_size;
    hsize_t     vnum_trial;
    hsize_t*    read_idx;
    
    ECoGMeta*   metadata;
    
    char        query_labels[QUERY_SIZE][LABEL_LEN];
    char        filename[HNAME_MAX];

    int         query_num, cio=1;
    int         proc_num, my_rank;
    
    double*     ECoGData_data;

    // Start MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    strcpy(filename, DATAFILE);

    // check arguments
    query_num = parse_argv(argc, argv, filename, query_labels, &cio);

    // Only root process open and read all index datasets
    if (my_rank == ROOTPROC) { 

        // Allocate matadata storage 
        metadata = (ECoGMeta*)malloc(sizeof(ECoGMeta));
        status = root_get_metadata(filename, metadata);
 
        // Allocate memory for actual read index
        read_idx = (hsize_t*)malloc(metadata->EIndx_ndims * sizeof(hsize_t));

        total_size = 0;
        // Find how much data we need to read
        for (i = 0; i < metadata->EIndx_dim[0]; i++) {
            for (j = 0; j < query_num; j++) {
                //printf("%s - %s\n", metadata->ELbls_data[ (int)(metadata->EIndx_data[i]) - 1 ], query_labels[j]);
                if ( strcmp( metadata->ELbls_data[ (int)(metadata->EIndx_data[i]) - 1 ], query_labels[j]) == 0 ) {
                    read_idx[total_size++] = i;
                }
            }
        }

        vnum_trial = metadata->vnum_trial;

        printf("Total trials to be read %llu\n", total_size);
    }
    
    // total_size is the number of data blocks(each with 301*256 elements) that we need to read
    // read_idx has all the offsets

    MPI_Bcast(&total_size, 1, MPI_LONG_LONG_INT, ROOTPROC, MPI_COMM_WORLD);
    
    // Allocate memory for actual read index, all same size for now
    if(my_rank != ROOTPROC)
        read_idx = (hsize_t*)malloc(total_size * sizeof(hsize_t));

    MPI_Bcast(read_idx, total_size, MPI_LONG_LONG_INT, ROOTPROC, MPI_COMM_WORLD);
    
    MPI_Bcast(&vnum_trial, 1, MPI_LONG_LONG_INT, ROOTPROC, MPI_COMM_WORLD);

        
    
    // All processes open the data file.
    plist_id        = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);

    file_id         = H5Fopen(filename, H5F_ACC_RDWR, plist_id);
    ECoGData_id     = H5Dopen(file_id, "/Data/ECoG", H5P_DEFAULT);
    ECoGData_space  = H5Dget_space(ECoGData_id);

    H5Sget_simple_extent_dims(ECoGData_space, ECoGData_dim, NULL);

    // All starts with 0 of the second dimension.    
    my_offset[1]       = 0;

    // Each reads vnum_trial(301) vectors.
    my_count[0]        = vnum_trial + 1;
    my_count[1]        = ECoGData_dim[1];

    // Total elements of each trial
    hsize_t elem_trial = my_count[1]*my_count[0];

    if (my_rank != proc_num-1) {
        my_trials = total_size/proc_num;
    }
    else {
        my_trials = total_size/proc_num + total_size%proc_num;
    }

    // Used for union selection
    hsize_t my_total_count[MAXDIM];
    my_total_count[0]  = my_count[0] * my_trials; 
    my_total_count[1]  = my_count[1]; 


    //ECoGData_memspace  = H5Screate_simple(2, my_count, NULL);
    ECoGData_memspace  = H5Screate_simple(2, my_total_count, NULL);

    // Actual data storage space
    ECoGData_data      = (double*)malloc(my_trials*elem_trial*sizeof(double));

    MPI_Barrier(MPI_COMM_WORLD);

    for (i = 0; i < my_trials; i++) {
        //printf("%d - my_i:%d\n", my_rank, (total_size/proc_num)*my_rank+i );
        my_offset[0] = read_idx[(total_size/proc_num)*my_rank+i] * my_count[0]; 

        if (i == 0)
            status   = H5Sselect_hyperslab(ECoGData_space, H5S_SELECT_SET, my_offset, NULL, my_count, NULL);
        else
            status   = H5Sselect_hyperslab(ECoGData_space, H5S_SELECT_OR, my_offset, NULL, my_count, NULL);
        //status = H5Sselect_hyperslab(ECoGData_space, H5S_SELECT_SET, my_offset, NULL, my_count, NULL);
        //printf("%d - Offsets: [0]:%llu [1]:%llu\n", my_rank, my_offset[0],my_offset[1]);
        //status = H5Dread(ECoGData_id, H5T_IEEE_F64LE, ECoGData_memspace, ECoGData_space, H5P_DEFAULT, ECoGData_data+elem_trial*i);
        //printf("%d@%d - %f %f\n", my_rank, i, ECoGData_data[elem_trial*(i+1) - 1], ECoGData_data[elem_trial*(i+1)-2]);
    }
 
    status = H5Dread(ECoGData_id, H5T_IEEE_F64LE, ECoGData_memspace, ECoGData_space, H5P_DEFAULT, ECoGData_data);
    if (status < 0)
        printf("Error reading!\n");
    
    double my_sum = 0.0, all_sum;
    for (i = 0; i < my_trials*elem_trial; i++) {
        // There are -nan exist in the datafile, so need to check
        if ( !isnan(ECoGData_data[i]) ) 
            my_sum += ECoGData_data[i]; 
    }
    
    //printf("\n");
    //printf("my total count: %llu %llu\n", my_total_count[0], my_total_count[1]);
    MPI_Reduce(&my_sum, &all_sum, 1, MPI_DOUBLE, MPI_SUM, ROOTPROC, MPI_COMM_WORLD);
    if(my_rank == ROOTPROC) {
        printf("My sum = %f, All sum = %f\n", my_sum, all_sum);

        /* for (i = my_trials*elem_trial; i > my_trials*elem_trial - 257 ; i--) { */
        /*     if(i % 6== 0) */
        /*         printf("\n"); */
        /*     printf("%f\t", ECoGData_data[i]); */
        /* } */

        printf("Start cleaning up\n");
    }

    // Close everything. 
    status = H5Dclose(ECoGData_id);
    status = H5Fclose(file_id);

    if (my_rank == ROOTPROC) {
        free(metadata->ECoGIndx_data);
        free(metadata->EIndx_data);
        for (i = 0; i < metadata->ELbls_dim[0]; i++) {
            free(metadata->ELbls_data[i]);
        }
        free(metadata->ELbls_data);
        free(metadata);
    }

    free(read_idx);
    free(ECoGData_data);

    H5close();
    MPI_Finalize();
    return 0;

}

herr_t root_get_metadata(char* filename, ECoGMeta* metadata) 
{
    int     i;
    hid_t   file_id;
    hid_t   ECoGIndx_id, EIndx_id, ELbls_id, ELbls_memtype;
    hid_t   ECoGIndx_space, EIndx_space, ELbls_space;
    hsize_t ECoGIndx_size=1, EIndx_size=1, ELbls_size=1;
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

    // Calculate how much space we need
    for (i = 0; i < metadata->ECoGIndx_ndims; i++) 
        ECoGIndx_size *= metadata->ECoGIndx_dim[i];

    for (i = 0; i < metadata->EIndx_ndims; i++) 
        EIndx_size    *= metadata->ECoGIndx_dim[i];

    // Allocate memory
    metadata->ECoGIndx_data = (double*)malloc(ECoGIndx_size*sizeof(double));
    metadata->EIndx_data    = (double*)malloc(EIndx_size*sizeof(double));
    
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

    metadata->vnum_trial = metadata->ECoGIndx_data[1] - metadata->ECoGIndx_data[0];

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

    return status;
}

int parse_argv(int argc, char* argv[], char* filename, char query_labels[QUERY_SIZE][LABEL_LEN], int* cio)
{
    int i = 0;
    char c;
    while ((c = getopt (argc, argv, "f:q:c:")) != -1) {
        switch (c) 
        {   
        case 'f':
            strcpy(filename, optarg);
            break;
        case 'q':
            if (i >= QUERY_SIZE ) {
                printf("Increase QUERY_SIZE in .h file\nExit\n");
                exit(-1);
            }
            strcpy(query_labels[i++], optarg);
            break;
        case 'c':
            *cio = atoi(optarg);
            break;
        default:
            printf("Error option [%s]\n", optarg);
            printf("%s -f filename -c collective_io -q query_label\n", argv[0]);
            exit(-1);
        }   
    }  
    return i;
}

int test(ECoGMeta* metadata)
{
    printf("Event_ECoGIndx: [0]:%.0f [1]:%.0f\n", metadata->ECoGIndx_data[0], metadata->ECoGIndx_data[1]);
    printf("Event_ECoGIndx: [3510]:%.0f [3511]:%.0f\n", metadata->ECoGIndx_data[3510], metadata->ECoGIndx_data[3511]);
    printf("EIndx_data: [0]:%.0f [1755]:%.0f\n", metadata->EIndx_data[0], metadata->EIndx_data[1755]);
    printf("Labels: [0]:%s [50]:%s\n", metadata->ELbls_data[0], metadata->ELbls_data[50]);

    return 0;
}
