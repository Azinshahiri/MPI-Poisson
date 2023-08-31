//--------------------
//  Header files
//--------------------
#include <mpi.h>
#include <math.h>
#include <iostream>
#include <fstream>
using namespace std;




//--------------------
//  Lattice size
//--------------------
const int N = 20;
int nx, ny;
const double h = 1.0 / (N+1);
const double h2 = h * h;
const double dim = 2;
const double degree = 2 * dim;
//--------------------
//  Iterations
//--------------------
const int iterations = 400;





//----------------------------------------
//      MAIN function:    u_global
//----------------------------------------
int main()
{

    
    //==================================================
    //  Preparation:    Exact solution + inhomogeneity
    //==================================================
    
    
    //--------------------
    //  Exact solution
    //--------------------
    double u_exact[ N + 2 ][ N + 2 ];
    for ( nx = 0; nx <= N + 1; nx++ )
    for ( ny = 0; ny <= N + 1; ny++ )
    {
        u_exact[nx][ny] = pow(h * nx,2) + pow(h * ny,2);
    }
    //------------------------------
    //  Alternative solutions
    //------------------------------
    // u_exact[nx][ny] = pow(h * nx,2) + pow(h * ny,2);
    // u_exact[nx][ny] = sin( 2*M_PI*nx / (N+1) ) * sin( 2*M_PI*ny / (N+1) );
    
    
    //--------------------
    //  Inhomogeneity
    //--------------------
    double f[N+2][N+2];
    for (nx=1; nx<=N; nx++)
    for (ny=1; ny<=N; ny++)
    {
        f[nx][ny] = 4;
    }
    //------------------------------
    //  Alternative inhomogeneities
    //------------------------------
    // f[nx][ny] = 4;
    // f[nx][ny] = - 8*( M_PI * M_PI ) * sin( 2*M_PI*nx / (N+1) ) * sin( 2*M_PI*ny / (N+1) );
    
    //==========================================================================================
    
    

    

    //--------------------
    // Initialize MPI
    //--------------------
    MPI_Init(NULL,NULL);
    const int MPI_root = 0;
    const int MPI_tag = 0;

    //--------------------
    //  MPI size and rank
    //--------------------
    int rank, processes;
    MPI_Comm_size(MPI_COMM_WORLD, &processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    //--------------------
    //  Interval length
    //--------------------
    const int partitions = 4;
    const int N_local = N / partitions;
    //------------------------------------------------------------
    //  Issue:  (parititions = const) = (processes = var)
    //
    //  =>  (array = fixed length) ≠ (pointer = dynamic alloc)  !!
    //------------------------------------------------------------




    //----------------------------------------
    //  Global array:   Dirichlet boundary
    //----------------------------------------
    double u_global[ N + 2 ][ N + 2 ] = {};
    if ( rank == MPI_root )
    {
        cout << "Global array: setting Dirichlet boundary" << endl;
        for ( ny = 0; ny <= N+1; ny++ )
        {
            u_global[0][ny]   = u_exact[0][ny];
            u_global[N+1][ny] = u_exact[N+1][ny];
        }
        for ( nx = 0; nx <= N+1; nx++ )
        {
            u_global[nx][0]   = u_exact[nx][0];
            u_global[nx][N+1] = u_exact[nx][N+1];
        }
    }
    //========================================




    //------------------------------------------------------------
    //  Local arrays:   [ 0 -- 1 -- N_local -- N_local + 1 ]
    //------------------------------------------------------------
    double u_local[ N_local + 2 ][ N + 2 ];


    int batchsize = ( N_local + 2 ) * ( N + 2 );
    int overlap = ( N + 2 );


    int sendcount[processes];
    int displacement[processes];
    for (int r = 0; r < processes; r++)
    {
        displacement[r] = r * (batchsize - 2 * overlap);
        sendcount[r] = batchsize;
    }
    MPI_Scatterv(u_global, sendcount, displacement, MPI_DOUBLE, u_local, batchsize , MPI_DOUBLE, MPI_root, MPI_COMM_WORLD);





    //============================================================
    //  Used in class:       Gauss Seidel
    //============================================================


    //----------------------------------------
    //  Iterations:     Gauss Seidell
    //----------------------------------------
    cout << "Starting Gaus Seidel:\t" << rank << endl;
    for (int iteration = 1; iteration <= iterations; iteration++)
    {

        //----------------------------------------
        //  Gauss Seidel:   [ 1 -- N_local ]
        //----------------------------------------
        for ( ny = 1; ny <= N; ny++ )
        for ( nx = 1; nx <= N_local; nx++ )
        {
            u_local[nx][ny]
            = ( - h2 * f[ nx + rank * N_local ][ny] + u_local[ nx-1 ][ny] + u_local[ nx+1 ][ny] + u_local[nx][ ny+1 ] + u_local[nx][ ny-1 ] ) /degree;
        }
        //----------------------------------------
        //  Boundary values:    send to right
        //----------------------------------------
        if ( rank < processes - 1 )
        {
            MPI_Send(&u_local[N_local], overlap, MPI_DOUBLE, rank + 1, MPI_tag, MPI_COMM_WORLD);
        }
        //----------------------------------------
        //  Boundary values:    send to left
        //----------------------------------------
        if ( rank > 0 )
        {
            MPI_Send(&u_local[1], overlap, MPI_DOUBLE, rank - 1, MPI_tag, MPI_COMM_WORLD);
        }
        //----------------------------------------
        //  Boundary values:    receive from right
        //----------------------------------------
        if ( rank < processes - 1 )
        {
            MPI_Recv(&u_local[N_local + 1], overlap, MPI_DOUBLE, rank + 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        //----------------------------------------
        //  Boundary values:    receive from left
        //----------------------------------------
        if ( rank > 0 )
        {
            MPI_Recv(&u_local[0], overlap, MPI_DOUBLE, rank - 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }
    //==============================






    //--------------------
    //  Synchronize MPI
    //--------------------
    MPI_Barrier(MPI_COMM_WORLD);




    //----------------------------------------
    //  Global array:   gather values
    //----------------------------------------
    if (rank==MPI_root)

    cout << "Global array: gathering local arrays" << endl;
    MPI_Gatherv(u_local, batchsize, MPI_DOUBLE, u_global, sendcount, displacement, MPI_DOUBLE, MPI_root, MPI_COMM_WORLD);



    //------------------------------
    //  Global array:   after gather
    //------------------------------
    if (rank == MPI_root)
    {
        cout << "Global array:" << endl;
        for ( ny=0; ny <= N + 1; ny++ )
        {
            for ( nx=0; nx <= N + 1; nx++ )
            {
                cout << u_global[nx][ny] << "\t";
            }
            cout << endl;
        }
    }
    //==============================



    //------------------------------------------------------------
    //  Signed difference:      u_diff = u_exact - u_global
    //------------------------------------------------------------
    double u_diff[ N + 2 ][ N + 2 ];
    if (rank == MPI_root)
    {
        cout << "Signed difference:" << endl;
        for ( ny=0; ny <= N + 1; ny++)
        {
            for ( nx=0; nx <= N + 1; nx++)
            {
                u_diff[nx][ny] = u_exact[nx][ny] - u_global[nx][ny];
                cout << u_diff[nx][ny] << "\t";
            }
            cout << endl;
        }
    }
    //==============================



    //------------------------------------------------------------
    //  Integration:        int \| u_exact - u_global \|^2 dx * dy
    //  Normalisation:      int     \| u_exact \|^2        dx * dy
    //------------------------------------------------------------
    if (rank == MPI_root)
    {
        const double differential = pow(h,dim);
        double integral = 0, normalisation = 0;
        for ( ny=0; ny <= N + 1; ny++)
        for ( nx=0; nx <= N + 1; nx++)
        {
            integral      += pow( u_diff[nx][ny],  2) * differential;
            normalisation += pow( u_exact[nx][ny], 2) * differential;
        }
        cout << "Normalised error:\t" << integral/normalisation << endl;
    }
    //==============================





    //------------------------------
    //  Saving results  >  text file
    //------------------------------
    if ( rank == MPI_root )
    {

        //--------------------
        //  Opening file
        //--------------------
        ofstream approx_file, diff_file;

        cout << "OS > opening files" << endl;
        approx_file.open("/Users/azinshahiri/Desktop/mpi-projects/mpi-poisson/mpi-poisson/approximation.txt");
        diff_file.open("/Users/azinshahiri/Desktop/mpi-projects/mpi-poisson/mpi-poisson/difference.txt");



        //--------------------
        //  Writing results
        //--------------------
        if ( approx_file.is_open() && diff_file.is_open() )
        {
            cout << "OS > writing results" << endl;
            for ( ny=0; ny <= N + 1; ny++ )
            {
                for ( nx=0; nx <= N + 1; nx++ )
                {
                    approx_file << u_global[nx][ny] << "\t";
                    diff_file << u_diff[nx][ny] << "\t";
                }
                approx_file << endl;
                diff_file << endl;
            }
        }
        else cout << "OS > opening failed" << endl;


        //------------------------------
        //  Closing files   (for savety)
        //------------------------------
        approx_file.close();
        diff_file.close();
        cout << "OS > files closed" << endl;
    }
    //==============================




    
    
    


    //--------------------
    //  Finalize MPI
    //--------------------
    MPI_Finalize();



}
//--------------------
//  END-OF:     main
//--------------------
