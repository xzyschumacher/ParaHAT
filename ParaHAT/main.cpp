/*
	Description: first parameter: lenght of seed 
	second parameter: choice of producing hash table
	third paramter: number of reference area
	forth paramter: whether writing cigar
	fifth: the path of comparision file;
*/
#include "aligner.h"

#include <time.h>
#include <stdio.h>
#include <iostream>
#include <cstdlib>
#include <sys/times.h>
#include"mpi.h"
using namespace std;

#define PRINT_LOG 

const std::string getCurrentDateTime() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
    // for more information about date/time format
    strftime(buf, sizeof(buf), "[INFO] %Y-%m-%dT%X", &tstruct);

    return buf;
}

int main(int argc, char *argv[])
{
	MPI_Init(&argc, &argv);

	int my_rank, comm_sz;

	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

	//clock_t start, finish;
	//double duration;
	long BeginTime = 0;
    long EndTime = 0;
	long tck = 0;
	struct tms tTMS;

	opts *opt = new opts;
	Form fm(opt);
	if (fm.opt_parse(argc,argv,opt)!=1)
		exit(1);
	
#ifdef PRINT_LOG
	fprintf(stderr,"%s [%s] started from my_rank = %d\n",getCurrentDateTime().c_str(), PACKAGE_NAME, my_rank);
#endif

    //start = clock();
	BeginTime = times(&tTMS);
	
	Aligner alig(opt);
	alig.Runtask();	

	if ( NULL != opt ) delete opt;

	//finish = clock();
	//duration = (double)(finish - start) / CLOCKS_PER_SEC;
	EndTime = times(&tTMS);
	tck = sysconf(_SC_CLK_TCK);

#ifdef PRINT_LOG
	fprintf(stderr,"%s [%s] ended from my_rank = %d running time = %f seconds\n",getCurrentDateTime().c_str(), PACKAGE_NAME, my_rank, (EndTime - BeginTime) / (double)tck);
#endif	

	MPI_Finalize();
	return 0;

}

