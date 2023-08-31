//--------------------
//  Header files
//--------------------
#include <mpi.h>
#include <iostream>
using namespace std;




//------------------------------
//  Exact solution:  u_exact
//------------------------------
//double* uexact(int N)
//{
//    double h, x;
//    double *uex;
//    int i;
//
//    h = 1.0 / N;
//    uex = calloc((N + 1), sizeof(double));
//
//    for (i = 0; i <= N; i++)
//    {
//        x = i * h;
//        uex[i] = 0.5 * x * (x - 1.0);
//    }
//    return uex;
//}


//--------------------------------------------------
//  Euclidean distance:     \|u_exact - u_approx\|
//--------------------------------------------------
//double euclidean_distance(int N, double *u_approx, double *uexact) {
//    double diff, error_sum;
//    int i;
//
//    error_sum = 0.0;
//    for (i = 0; i <= N; i++)
//    {
//        diff = uexact[i] - u_approx[i];
//        error_sum += diff * diff;
//    }
//    return sqrt(error_sum);
//}



//----------------------------------------
//      MAIN function:    u_approx
//----------------------------------------
int main()
{


    //------------------------------
    //  Local approximation
    //------------------------------
    double *u_local;
    //------------------------------
    //  Iteration counters
    //------------------------------
    int i, j;
    //--------------------
    //  Lattice size
    //--------------------
    int N = 32;
    int h = 1.0 / N;
    int h2 = h * h;
    //--------------------
    //  Boundary values
    //--------------------
    double a = 0.0;
    double b = 0.0;
    //--------------------
    //  Iterations
    //--------------------
    int outer = 200;
    int inner = 6;
    


    //--------------------
    // Initialize MPI
    //--------------------
    MPI_Init(NULL,NULL);
    //--------------------
    //  MPI size and rank
    //--------------------
    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    //--------------------
    //  Domain intervals
    //--------------------
    int N_local = N / size;



    //------------------------------
    //  Global array:   u_approx
    //------------------------------
    if (rank == 0)
    {
        u_approx = calloc((N + 1), sizeof(double));
        u_approx[0] = a;
        u_approx[N] = b;
    }
    //----------------------------------------------------------------------
    //  Local arrays:   u_local:    [ 0 -- 1 -- N_local -- N_local + 1 ]
    //----------------------------------------------------------------------
    u_local = (double *) calloc(N_local + 2, sizeof(double));
    
    if (rank == 0)
    u_local[0] = a;

    if (rank == size - 1)
    u_local[N_local + 1] = b;
    

    
    //----------------------------------------
    //  Gauss Seidel:   Outer iterations
    //----------------------------------------
    for (int k = 0; k < outer; k++)
    {

        //--------------------------------------------------
        //  Inner iterations
        //--------------------------------------------------
        //  Gauss Seidel > inner points:    [ 1 -- N_local ]
        //--------------------------------------------------
        for ( j = 0; j < inner; j++ )
        for ( i = 1; i <= N_local; i++ )
        {
            u_local[i] = (-h2 + u_local[i - 1] + u_local[i + 1]) / 2.0;
        }
        //------------------------------
        //  Update:     Left boundary
        //------------------------------
        if ( rank < size - 1 )
        {
            MPI_Send(&u_local[N_local], 1, MPI_DOUBLE, rank + 1, MPI_ANY_TAG, MPI_COMM_WORLD);
            MPI_Recv(&u_local[N_local + 1], 1, MPI_DOUBLE, rank + 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        //------------------------------
        //  Update:     Right boundary
        //------------------------------
        if ( rank > 0 )
        {
            MPI_Send(&u_local[1], 1, MPI_DOUBLE, rank - 1, MPI_ANY_TAG, MPI_COMM_WORLD);
            MPI_Recv(&u_local[0], 1, MPI_DOUBLE, rank - 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

    }

    
    //--------------------
    //  Domain assembly
    //--------------------
    int root = 0;
    int* u_approx;
    u_approx = (double *) calloc(N+2, sizeof(double));

    MPI_Gather(&u_local[1], N_local, MPI_INT, &u_approx[1], N_local, MPI_INT, root, MPI_COMM_WORLD);
    u_approx[0]=a;
    u_approx[N+2]=B;



    //--------------------
    //  Exact solution
    //--------------------
    int* u_exact;
    u_exact = uexact(N);



    //----------------------------------------
    //      Comparison:  exact / approx
    //----------------------------------------
    cout << "Comparison: exact / approx" << str::endl;
    for (i = 0; i <= N; i++)
    {
        printf("%d\t %f\t %f\n", i, u_approx[i], u_exact[i] );
    }
    //----------------------------------------
    //      Distance:  exact / approx
    //----------------------------------------
    double error = euclidean_distance(N, u_approx, u_exact);

    printf("Distance:\t %f in iteration: %d\n", error, outer * inner);
    printf("value of u in midpoint: %f\n", u_approx[N / 2]);
    printf("exact value of u in midpoint: %f\n", u_exact[N / 2]);
    printf("number of intervals: %d for process:%d\n", N, rank);


    //--------------------
    //  Finalize MPI
    //--------------------
    MPI_Finalize();


    //--------------------
    //  Free up memory
    //--------------------
    free(u_approx);
    free(u_local);
    free(u_temp);
    free(u_exact);
    return 0;
}
