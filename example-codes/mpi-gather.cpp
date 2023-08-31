/*===============
    IO stream
 ===============*/
#include <iostream>
using namespace std;



/*====================
    Random numbers
 ====================*/
#include <stdlib.h>     // rand, srand



/*====================
    Print functions
 ====================*/
#define print(x) cout << x << endl
void outcome(int rank, int value)
{
    printf("rank: %d\t value: %d\n",rank,value);
}



/*====================
    Pausing process
 ====================*/
#include <chrono> // seconds
#include <thread> // sleep
using namespace std::chrono;
using namespace std::this_thread;



/*===============
    Open MPI
 ===============*/
#include <mpi.h>



/*===============
    Fill array
 ===============*/
void fill(int array[], int length)
{
    for (int i=0; i<length; i++)
    array[i] = rand()%20 ;
}



/*====================
    Compute average
 ====================*/
int average(int array[], int length)
{
    int sum=0;
    
    for (int i=0; i<length; i++)
    sum += array[i];
    
    return sum/length;
}



/*====================
    Main function
 ====================*/
int main()
{

    
    /*====================
        Initialize MPI
     ====================*/
    MPI_Init(NULL, NULL);
    /*=========================
        MPI: size and rank
     =========================*/
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    
    //------------------------------
    //  Subarrays initialise
    //------------------------------
    int sublength = 3+2;
    int* locarray;
    locarray = (int *) calloc(sublength, sizeof(int));
    //--------------------
    //  Setting subarrays
    //--------------------
    for (int i=0; i < sublength; i++)
    {
        locarray[i] = i + rank * sublength;
    }
    //--------------------
    //  Printing subarrays
    //--------------------
    printf("Subarray:\t %d\t>\t",rank);
    for (int i=0; i < sublength; i++)
    {
        cout << locarray[i] << "\t";
    }
    cout << endl;
    
    
    
    //------------------------------
    //  Mainarray:  offsets
    //------------------------------
    int root = 0;
    int suboffset = 1;
    int submany = sublength - suboffset-1;
    int mainoffset = 1;
    int mainextra = 1;
    int min_length = size * submany + mainoffset + mainextra;
    
    //------------------------------
    //  Mainarray   >  MPI gather
    //------------------------------
    cout << submany << "\t" << min_length << endl;
    int* mainarray;
    mainarray = (int *) calloc(min_length, sizeof(int));
    MPI_Gather(&locarray[suboffset], submany, MPI_INT, &mainarray[mainoffset], submany, MPI_INT, root, MPI_COMM_WORLD);
    
    //--------------------
    //  Display mainarray
    //--------------------
    if (rank == root)
    {
        cout << "Mainarray:" << endl;
        for (int j=0; j < min_length; j++)
        {
            cout << mainarray[j] << "\t";
        }
        cout << endl;
    }
    
    
    
    
    /*=========================
        Array + Batch size
     =========================*/
    int batch = 3;
    int length = batch*size;
    int array[length];
    int subarray[batch];
    /*===============
        Fill array
     ===============*/
    if (rank==0)
    {
        fill(array,length);
        
        cout << "Random numbers:" << endl;
        for (int i=0; i<length; i++)
        cout << array[i] << "\t";
        cout << endl;
        
        average(array, length);
    }
    /*=========================
        END-OF: preparation
     =========================*/
    
    
    
    /*====================
        Scatter array
     ====================*/
    MPI_Scatter(&array, batch, MPI_INT, &subarray, batch, MPI_INT, 0, MPI_COMM_WORLD);
    
    
    
    /*====================
        Batch averages
     ====================*/
    int batch_average = average(subarray, batch);
    outcome(rank,batch_average);
    
    
    

    /*====================
        Gather averages
     ====================*/
    int averages[size];
    MPI_Gather(&batch_average, 1, MPI_INT, &averages, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    
    
   /*====================
        Total average
    ====================*/
    if (rank==0)
    {
        int total_average = average(averages, size);
        cout << "Total average:\t" << total_average << endl;
    }
    /*====================
        Finalize MPI
     ====================*/
    MPI_Finalize();
}
