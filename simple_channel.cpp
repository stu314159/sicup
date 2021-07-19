#include "palabos2D.h"
#include "palabos2D.hh"

using namespace plb;
using namespace plb::descriptors;
using namespace std;

typedef double T;
#define DESCRIPTOR descriptors::D3Q19Descriptor

int main(int argc, char* argv[])
{
	plbInit(&argc,&argv);
	plint numCores = global::mpi().getSize();
	pcout << "Number of MPI processes: " << numCores << std::endl;

	plint N;
	try {
		global::argv(1).read(N);
	}
	catch(...)
	{
		pcout << "Wrong parameters. The syntax is " << std::endl;
		pcout << argv[0] << " N" << std::endl;
		pcout << "where N is the resolution of the mesh " << std::endl;
		exit(1);
	}

	// create incompressible flow parameters object
	IncomprFlowParam<T> parameters(
			(T) 1e-2, // uMax
			(T) 100., // Re
			N,        // N (read from input
			1.,      // lx
			1.,      // ly
			5.       // lz
	);

	return 0;
}
