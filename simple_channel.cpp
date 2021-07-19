#include "palabos2D.h"
#include "palabos2D.hh"
#include "palabos3D.h"
#include "palabos3D.hh"

#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>

// essentially a copy/adaptation of the channelWithCylinder3d example
using namespace plb;
using namespace plb::descriptors;
using namespace std;

typedef double T;
#define DESCRIPTOR descriptors::D3Q19Descriptor

#define NMAX 150

const T pi = (T)4.*std::atan((T)1.);


struct Param // parameters for turbulence model and "sponge zone"
{
	plint numOutletSpongeCells; // reference for Sponge cells?
	int outletSpongeZoneType;
	T targetSpongeCSmago;
};

Param param;

// object that will determine if a lattice point is inside or outside of the cylinderical shape.
// it pubically inherits the functional class so it will overload operator()
template<typename T>
class CylinderShapeDomain3D : public plb::DomainFunctional3D 
{
	public: 
		CylinderShapeDomain3D(plb::plint cx_, plb::plint cz_, plb::plint radius)
			: cx(cx_), cz(cz_), radiusSqr(plb::util::sqr(radius)) { }
		virtual bool operator()(plb::plint iX, plb::plint iY, plb::plint iZ) const {
			return plb::util::sqr(iX-cx) + plb::util::sqr(iZ-cz) <= radiusSqr;
		}
		virtual CylinderShapeDomain3D<T>* clone() const {
			return new CylinderShapeDomain3D<T>(*this);
		}
	private:
		plb::plint cx;
		plb::plint cz;
		plb::plint radiusSqr;
};

// a pressure boundary condition.  
static T poiseuillePressure(IncomprFlowParam<T> const &parameters, plint maxN)
{
	const T a =  parameters.getNx() - 1;
	const T b =  parameters.getNy() - 1;
	const T nu = parameters.getLatticeNu();
	const T uMax = parameters.getLatticeU();

	T sum = T(); // I guess this initializes to 0?
	for (plint iN = 0; iN < maxN; iN += 2)
	{
		T twoNplusOne = (T) 2*(T)iN + (T)1;
		sum += ((T)1/(std::pow(twoNplusOne,(T)3)*std::cosh(twoNplusOne*pi*b/((T)2*a))));
	}
	for (plint iN = 1; iN < maxN; iN +=2)
	{
		T twoNplusOne = (T)2*(T)iN + (T)1;
		sum -= ((T) 1/(std::pow(twoNplusOne,(T)3)*std::cosh(twoNplusOne*pi*b/((T)2*a))));
	}
	T alpha = -(T)8 * uMax*pi*pi*pi/(a*a*(pi*pi*pi - (T)32*sum)); //alpha = -dp/dz/mu
	T deltaP = -(alpha*nu);

	return deltaP;
}

T poiseuilleVelocity(plint iX, plint iY, IncomprFlowParam<T> const& parameters, plint maxN)
{
	const T a = parameters.getNx() - 1;
	const T b = parameters.getNy() - 1;

	const T x = (T)iX - a / (T)2;
	const T y = (T)iY - b / (T)2;

	const T alpha = -poiseuillePressure(parameters, maxN)/parameters.getLatticeNu();

	T sum = T();

	for (plint iN = 0; iN < maxN; iN += 2)
	{
		T twoNplusOne = (T)2*(T)iN + (T)1;
		sum += (std::cos(twoNplusOne*pi*x/a)*std::cosh(twoNplusOne*pi*y/a)
				/ (std::pow(twoNplusOne,(T)3)*std::cosh(twoNplusOne*pi*b/((T)2*a))));

	}
	for (plint iN = 1; iN < maxN; iN += 2)
	{
		T twoNplusOne = (T)2*(T)iN+(T)1;
		sum -= (std::cos(twoNplusOne*pi*x/a)*std::cosh(twoNplusOne*pi*y/a)
				/(std::pow(twoNplusOne,(T)3)*std::cosh(twoNplusOne*pi*b/((T)2*a))));
	}
	sum *= ((T)4 * alpha * a * a/std::pow(pi,(T)3));
	sum += (alpha / (T)2 * (x * x - a*a / (T)4));

	return sum;
}

template <typename T>
class SquarePoiseuilleDensityAndVelocity {
	public: 
		SquarePoiseuilleDensityAndVelocity(IncomprFlowParam<T> const& parameters_, plint maxN_)
			: parameters(parameters_), maxN(maxN_)
		{ }
		void operator()(plint iX, plint iY, plint iZ, T & rho, Array<T,3>& u) const {
			rho = (T)1;  // ^^ where is the Array template defined?
			u[0] = T();
			u[1] = T();
			u[2] = poiseuilleVelocity(iX,iY,parameters,maxN);
		}
	private:
		IncomprFlowParam<T> parameters;
		plint maxN;
};

template <typename T>
class SquarePoiseuilleVelocity {
	public:
		SquarePoiseuilleVelocity(IncomprFlowParam<T> const& parameters_, plint maxN_)
			: parameters(parameters_), maxN(maxN_)
		{}
		void operator()(plint iX, plint iY, plint iZ, Array<T,3>& u) const {
			u[0] = T();
			u[1] = T();
			u[2] = poiseuilleVelocity(iX,iY,parameters,maxN);
		}
	private:
		IncomprFlowParam<T> parameters;
		plint maxN;
};


void simulationSetup( MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
		IncomprFlowParam<T> const& parameters,
		OnLatticeBoundaryCondition3D<T,DESCRIPTOR>& boundaryCondition,
		Array<plint,3> & forceIds)
{
	// no periodic boundaries
	lattice.periodicity().toggleAll(false);

	const plint nx = parameters.getNx();
	const plint ny = parameters.getNy();
	const plint nz = parameters.getNz();

	// Flow direction: z
	// // Top/Bottom direction: y
	
	Box3D top = Box3D(0,nx-1,ny-1,ny-1,0,nz-1); // xyz coordinates of two corners of box
	Box3D bottom = Box3D(0,nx-1,0,0,0,nz-1); // xyz coordinates of two corners of lower surface box

	Box3D inlet = Box3D(0,nx-1,0,ny-1,0,0);
	Box3D outlet = Box3D(1,nx-2,1,ny-2,nz-1,nz-1); // note choice of x and y lattice point

	Box3D front = Box3D(0,0,0,ny-1,0,nz-1); // example calls this "right" 
	Box3D back = Box3D(nx-1,nx-1,0,ny-1,0,nz-1); // example calls this "left"

	boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice,inlet);
	boundaryCondition.addVelocityBoundary2P(outlet,lattice,boundary::neumann);

	boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice,top);
	boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice,bottom);

	boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice,back);
	boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice,front);

	setBoundaryVelocity(lattice,inlet,SquarePoiseuilleVelocity<T>(parameters,NMAX));
	setBoundaryVelocity(lattice,top,Array<T,3>((T)0.0, (T)0.0, (T)0.0));
	setBoundaryVelocity(lattice,bottom,Array<T,3>((T)0.0, (T)0.0, (T)0.0));
	setBoundaryVelocity(lattice,back,Array<T,3>((T)0.0, (T)0.0, (T)0.0));
	setBoundaryVelocity(lattice,front,Array<T,3>((T)0.0, (T)0.0, (T)0.0));

	// sponge zone parameters
	if (param.numOutletSpongeCells > 0)
	{
		T bulkValue;
		Array<plint,6> numSpongeCells;

		if (param.outletSpongeZoneType == 0) {
			pcout << "Generating an outlet viscosity sponge zone." << std::endl;
			bulkValue = parameters.getOmega();
		}
		else if (param.outletSpongeZoneType == 1) {
			pcout << "Generating an outlet Smagorinsky sponge zone. " << std::endl;
			bulkValue = 0.14;
		}
		else {
			pcout << "Error: unknown type of sponge zone." << std::endl;
			exit(-1);
		}

		// Number of sponge zone lattice nodes at all the outer domain boundaries
		// 0 means boundary at x = 0
		// 1 means boundary at x = nx-1
		// 2 means boundary at y = 0
		// etc...
		numSpongeCells[0] = 0;
		numSpongeCells[1] = 0;
		numSpongeCells[2] = 0;
		numSpongeCells[3] = 0;
		numSpongeCells[4] = 0;
		numSpongeCells[5] = param.numOutletSpongeCells;

		std::vector<MultiBlock3D*> args; // I guess this is why we need <vector>
		args.push_back(&lattice);

		if (param.outletSpongeZoneType == 0) {
			applyProcessingFunctional(new ViscositySpongeZone3D<T,DESCRIPTOR>(
						nx,ny,nz,bulkValue,numSpongeCells),
					lattice.getBoundingBox(),args);
		}
		else {
			applyProcessingFunctional(new SmagorinskySpongeZone3D<T,DESCRIPTOR>(
						nx,ny,nz,bulkValue,param.targetSpongeCSmago,numSpongeCells),
					lattice.getBoundingBox(),args);
		}
	}

	initializeAtEquilibrium(lattice,lattice.getBoundingBox(),SquarePoiseuilleDensityAndVelocity<T>(parameters,NMAX));
	// add the obstacle: cylidner 3d
	plint cx = nx/2 + 2; // to put the cylinder off the center of the x-dim of the channel
	plint cz = nz/6.;
	plint radius = parameters.getResolution() / 2.; // getResolution gives the reference length; obviously the radius is half of that in lattice units
	lattice.toggleInternalStatistics(true); // what does this do?
	forceIds[0] = lattice.internalStatSubscription().subscribeSum();
	forceIds[1] = lattice.internalStatSubscription().subscribeSum();
	forceIds[2] = lattice.internalStatSubscription().subscribeSum();

	defineDynamics(lattice, lattice.getBoundingBox(),
			new CylinderShapeDomain3D<T>(cx,cz,radius),
			new plb::MomentumExchangeBounceBack<T,DESCRIPTOR>(forceIds));
	initializeMomentumExchange(lattice,lattice.getBoundingBox());

	lattice.initialize();
}

template<class BlockLatticeT>
void writeGifs(BlockLatticeT& lattice,
		IncomprFlowParam<T> const& parameters, plint iter)
{
	const plint imSize = 600;
	const plint nx = parameters.getNx();
	const plint ny = parameters.getNy();
	const plint nz = parameters.getNz();

	Box3D slice(0,nx-1,ny/2,ny/2,0,nz-1); // length of channel sliced in 2 along the y-axis
	ImageWriter<T> imageWriter("leeloo");

	imageWriter.writeScaledGif(createFileName("uNorm",iter,6),
			*computeVelocityNorm(lattice,slice),
			imSize,imSize);
}

template<class BlockLatticeT>
void writeVTK(BlockLatticeT& lattice, IncomprFlowParam<T> const& parameters, plint iter)
{
	T dx = parameters.getDeltaX();
	T dt = parameters.getDeltaT();
	ParallelVtkImageOutput3D<T> vtkOut(createFileName("vtk",iter,6),3,dx);
	vtkOut.writeData<3,float>(*computeVelocity(lattice),"velocity",dx/dt);
	vtkOut.writeData<float>(*computeVelocityNorm(lattice),"velocityNorm",dx/dt);
	vtkOut.writeData<3,float>(*computeVorticity(*computeVelocity(lattice)),"vorticity",1./dt);

}

int main(int argc, char* argv[])
{
	plbInit(&argc,&argv);
	plint numCores = global::mpi().getSize();
	pcout << "Number of MPI processes: " << numCores << std::endl;
	
	global::directories().setOutputDir("./tmp/");

	// read command line parameters
	if (argc != 10)
	{
		pcout << "Error: the parameters are wrong. \n";
		pcout << "Give: Re Ndivs cylD Lx_p Ly_p Lz_p outletSpongeZoneType numOutletSpongeCells targetSpongeCSmago \n";
		pcout << "Example: ./simple_channel 100 20 6 50 50 75 1 1 0.6 \n";
		exit(1);
	}

	T Re_ = atof(argv[1]);
	plint N_ = atoi(argv[2]);
	T D_ = atof(argv[3]);
	T W_ = atof(argv[4]);
	T h_ = atof(argv[5]);
	T L_ = atof(argv[6]);

	param.outletSpongeZoneType = atoi(argv[7]);
	param.numOutletSpongeCells = atoi(argv[8]);
	param.targetSpongeCSmago = atof(argv[9]);

	// create incompressible flow parameters object
	IncomprFlowParam<T> parameters( 
			(T) 0.1, // uMax - must be in lattice units
			Re_, // Re
			N_,        // N (read from input) - resolution of the reference length
			W_/D_,      // lx (must be in dimensionless units!!)
			h_/D_,      // ly
			L_/D_       // lz
	);
	const T vtkSave = (T) 1/(T) 100; // sets periodicity of VTK output
	const T maxT = (T) 200.0; // sets maximum time of the simulation
	
	pcout << "omega = " << parameters.getOmega() << std::endl;
	writeLogFile(parameters, "3D square Poiseuille with Cylinder as an obstacle");

	T omega = parameters.getOmega();

	MultiBlockLattice3D<T,DESCRIPTOR> lattice(
		parameters.getNx(), parameters.getNy(), parameters.getNz(),
		new BGKdynamics<T,DESCRIPTOR>(omega));

	OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* boundaryCondition = 
		createLocalBoundaryCondition3D<T,DESCRIPTOR>();

	pcout << "With current settings, number of time steps: " << parameters.nStep(maxT) << std::endl;

	plint nnodes = parameters.getNx()*parameters.getNy()*parameters.getNz();
	pcout << "Number of lattice points: " << nnodes << std::endl;
	
	Array<plint,3> forceIds;
	simulationSetup(lattice,parameters,*boundaryCondition,forceIds);

	// loop over main time iteration
	for(plint iT=0; iT<parameters.nStep(maxT);++iT)
	{
		if(iT%parameters.nStep(vtkSave)==0 && iT>0)
		{
			pcout << "step " << iT << "; t=" << iT*parameters.getDeltaT() << std::endl;
			//writeGifs(lattice,parameters,iT);
			writeVTK(lattice,parameters,iT);
		}

		// execute a time iteration
		lattice.collideAndStream();

	}

	delete boundaryCondition;
	return 0;
}
