#include "palabos3D.h"
#include "palabos3D.hh"
#include "headersJo.hh"
#include "headersJo.h"

#include <cstdlib>
#include <iostream>

using namespace plb;
using namespace std;

typedef double T;


#define NSDESCRIPTOR descriptors::ForcedD3Q19Descriptor
#define ADESCRIPTOR descriptors::AdvectionDiffusionD3Q7Descriptor
#define NSDYNAMICS GuoExternalForceBGKdynamics
/// Initialization of the volume fraction field.
template<typename T, template<typename NSU> class nsDescriptor, template<typename ADU> class adDescriptor>
struct IniVolFracProcessor3D : public BoxProcessingFunctional3D_S<T>
{
    IniVolFracProcessor3D(RayleighTaylorFlowParam<T,nsDescriptor,adDescriptor> parameters_)
        : parameters(parameters_)
    { }
    virtual void process(Box3D domain, ScalarField3D<T>& volfracField)
    {
        Dot3D absoluteOffset = volfracField.getLocation();


      T nz=parameters.getNz();
      T nx=parameters.getNx();
      T ny=parameters.getNy();
      T up=0.135;
      T low=0.25;

      T kx;
      T ky;
            for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
                for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                    for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                      plint absoluteX = absoluteOffset.x + iX;
                      plint absoluteY = absoluteOffset.y + iY;
                        plint absoluteZ = absoluteOffset.z + iZ;

                        	T VolFrac;
                          T rand_val=(double)rand()/RAND_MAX;
                          T amplitude=0.5;

                          kx=(2*3.1416)/(nx/1);
                          ky=(2*3.1416)/(ny/3);


    			if (absoluteZ==(nz-1)-(int)((up/(up+low))*nz)) {
            VolFrac=1-(rand_val*amplitude*(1+cos(absoluteX*kx)*cos(absoluteY*ky)));
          }
          else
    			if(absoluteZ<=(nz-1)-(int)((up/(up+low))*nz)-1)
    				VolFrac=0.0;
			else
			if(absoluteZ==(nz-1))
			VolFrac=0.0;
    	else VolFrac=1.;



                    volfracField.get(iX, iY, iZ) = VolFrac;

                }
            }
        }
    }

    virtual IniVolFracProcessor3D<T,nsDescriptor,adDescriptor>* clone() const
    {
        return new IniVolFracProcessor3D<T,nsDescriptor,adDescriptor>(*this);
    }

    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::staticVariables;
    }
    virtual BlockDomain::DomainT appliesTo() const {
        return BlockDomain::bulkAndEnvelope;
    }
private :
    RayleighTaylorFlowParam<T,nsDescriptor,adDescriptor> parameters;
};

/// Initialization of the density field.
template<typename T, template<typename NSU> class nsDescriptor, template<typename ADU> class adDescriptor>
struct IniDensityProcessor3D : public BoxProcessingFunctional3D_S<T>
{
    IniDensityProcessor3D(RayleighTaylorFlowParam<T,nsDescriptor,adDescriptor> parameters_)
        : parameters(parameters_)
    { }
    virtual void process(Box3D domain, ScalarField3D<T>& densityField)
    {
        Dot3D absoluteOffset = densityField.getLocation();


        T nz=parameters.getNz();

        T up=0.135;
        T low=0.25;

              for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
                  for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                      for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {

                          plint absoluteZ = absoluteOffset.z + iZ;

                          	T dens;


            if(absoluteZ==0)
              dens=0.;
            else
      			if(absoluteZ<=(nz-1)-(int)((up/(up+low))*nz))
      				dens= 1008.4;
      			else
            if(absoluteZ==(nz-1))
            dens=0.;
            else
              dens=1000.;



                    densityField.get(iX, iY, iZ) = dens;
                }
            }
        }
    }
    virtual IniDensityProcessor3D<T,nsDescriptor,adDescriptor>* clone() const
    {
        return new IniDensityProcessor3D<T,nsDescriptor,adDescriptor>(*this);
    }

    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::staticVariables;
    }
    virtual BlockDomain::DomainT appliesTo() const {
        return BlockDomain::bulkAndEnvelope;
    }
private :
    RayleighTaylorFlowParam<T,nsDescriptor,adDescriptor> parameters;
};

/*
T poiseuilleDensity(plint iX, plint iY, plint iZ, RayleighTaylorFlowParam<T,NSDESCRIPTOR,ADESCRIPTOR> const& parameters) {
    T nz = parameters.getNz();
    T nx = parameters.getNx();
    T ny = parameters.getNy();
    T dt = parameters.getDeltaT();
    T dx = parameters.getDeltaX();
    T rhophys = 1000;
    T gphys = 9.81;
    if(iZ==0 || iZ==nz-1) return 0.;
    else
    if(iX==0 || iX==nx-1) return 0.;
    else
    if(iY==0 || iY==ny-1) return 0.;
    else
    return 1+(NSDESCRIPTOR<T>::invCs2*((dt*dt*rhophys*gphys*(nz-1-iZ))/dx));
}

template<typename T>
class PoiseuilleDensity {
public:
    PoiseuilleDensity(RayleighTaylorFlowParam<T,NSDESCRIPTOR,ADESCRIPTOR> parameters_)
        : parameters(parameters_)
    { }
    T operator()(plint iX, plint iY, plint iZ) const {
        return poiseuilleDensity(iX, iY, iZ, parameters);
    }
private:
    RayleighTaylorFlowParam<T,NSDESCRIPTOR,ADESCRIPTOR> parameters;
};

template<typename T>
class PoiseuilleDensityAndZeroVelocity {
public:
    PoiseuilleDensityAndZeroVelocity(RayleighTaylorFlowParam<T,NSDESCRIPTOR,ADESCRIPTOR> parameters_)
        : parameters(parameters_)
    { }
    void operator()(plint iX, plint iY, plint iZ, T& rho, Array<T,3>& u) const {
        rho = poiseuilleDensity(iX, iY, iZ,parameters);
        u[0] = T();
        u[1] = T();
        u[2] = T();
    }
private:
    RayleighTaylorFlowParam<T,NSDESCRIPTOR,ADESCRIPTOR> parameters;
};
*/


void ExpSetup (
        MultiBlockLattice3D<T, NSDESCRIPTOR>& nsLattice,
        MultiScalarField3D<T>& volfracField,
        MultiScalarField3D<T>& densityField,
        RayleighTaylorFlowParam<T,NSDESCRIPTOR,ADESCRIPTOR> &parameters )
{

//  plint nx = parameters.getNx();
//  plint ny = parameters.getNy();
//  plint nz = parameters.getNz();




  initializeAtEquilibrium(nsLattice, nsLattice.getBoundingBox(), (T)1., Array<T,3>((T)0.,(T)0.,(T)0.) );



    applyProcessingFunctional(
            new IniVolFracProcessor3D<T,NSDESCRIPTOR,ADESCRIPTOR>(parameters),
            volfracField.getBoundingBox(), volfracField );

    applyProcessingFunctional(
            new IniDensityProcessor3D<T,NSDESCRIPTOR,ADESCRIPTOR>(parameters),
            densityField.getBoundingBox(), densityField );


    nsLattice.initialize();

}




void writeVTK(MultiBlockLattice3D<T,NSDESCRIPTOR>& nsLattice,
              MultiScalarField3D<T>& volfracField,
              MultiScalarField3D<T>& densityField,
              RayleighTaylorFlowParam<T,NSDESCRIPTOR,ADESCRIPTOR> const& parameters, plint iter)
{

    T dx = parameters.getDeltaX();
    T dt = parameters.getDeltaT();

    VtkImageOutput3D<T> vtkOut(createFileName("vtk", iter, 6), dx);
    //vtkOut.writeData<float>(*computeVelocityNorm(nsLattice), "NSvelocityNorm", dx/dt);
    vtkOut.writeData<3,float>(*computeVelocity(nsLattice), "NSvelocity", dx/dt);
    //vtkOut.writeData<3,float>(*computeVorticity(*computeVelocity(nsLattice)), "vorticity", 1./dt);
    // Temperature is the order-0 moment of the advection-diffusion model. It can
    //    therefore be computed with the function "computeDensity".
    vtkOut.writeData<float>(volfracField, "VolFrac", (T)1);

}

void writeGif(MultiBlockLattice3D<T,NSDESCRIPTOR>& nsLattice,
              MultiScalarField3D<T>& volfracField,
              MultiScalarField3D<T>& densityField,
	      MultiScalarField3D<T>& v_sed, int iT)
{
    const plint imSize = 600;
    const plint nx = nsLattice.getNx();
    const plint ny = nsLattice.getNy();
    const plint nz = nsLattice.getNz();
    Box3D slice((nx-1)/2, (nx-1)/2, 0, ny-1, 0, nz-1);
    ImageWriter<T> imageWriter("water.map");
    imageWriter.writeScaledGif(createFileName("u", iT, 6),
                               *computeVelocityNorm(nsLattice, slice),
                               imSize, imSize);

    // Temperature is the order-0 moment of the adveadction-diffusion model. It can
    //    therefore be computed with the function "computeDensity".
    imageWriter.writeScaledGif(createFileName("VolFrac", iT, 6),
                               *extractSubDomain(volfracField, slice),
                               imSize, imSize);

    imageWriter.writeScaledGif(createFileName("Density", iT, 6),
                                *extractSubDomain(densityField, slice),
                                imSize, imSize);

    imageWriter.writeScaledGif(createFileName("v_sed", iT, 6),
                                *extractSubDomain(v_sed, slice),
                                imSize, imSize);

    imageWriter.writeScaledGif( createFileName("Vorticity", iT, 6),
                                *computeNorm(*computeVorticity (
                                        *computeVelocity(nsLattice) ), slice ),
                                imSize, imSize );
}


int main(int argc, char *argv[])
{
    plbInit(&argc, &argv);
    global::timer("simTime").start();




    const T lx  = 0.075;
    const T ly  = 0.303;
    const T lz  = 0.385;
    const T uCar  = 0.008;
    const T uMax = 3.2e-4;
    const T Gr = 1e6;
    //const T Kappa = 1e-6;
    const T Ri = 0.10541;
    const T Di = 3e-9;

    const plint resolution = 500;

    global::directories().setOutputDir("./tmp/");



    RayleighTaylorFlowParam<T,NSDESCRIPTOR,ADESCRIPTOR> parameters (
            Ri,
            Gr,
            uMax,
            uCar,
            resolution,
            lx,
            ly,
            lz,
            Di );



    T rho0=1000;
    T rhoP=2550;
    T tunedg=9.81*parameters.getLatticeGravity();
    T rho0f=0.0;
    T omega;



    writeLogFile(parameters,"palabos.log");


    plint nx = parameters.getNx();
    plint ny = parameters.getNy();
    plint nz = parameters.getNz();


    T nsOmega = parameters.getSolventOmega();
    T convers = parameters.getDeltaT()/parameters.getDeltaX();

    Dynamics<T,NSDESCRIPTOR>* nsdynamics = new NSDYNAMICS<T,NSDESCRIPTOR>(nsOmega);
    MultiBlockLattice3D<T,NSDESCRIPTOR> nsLattice(nx,ny,nz, nsdynamics->clone());
    defineDynamics(nsLattice, nsLattice.getBoundingBox(), nsdynamics->clone());
    delete nsdynamics; nsdynamics = 0;



    Box3D initDomain(nsLattice.getBoundingBox());
    initDomain.x0++;
    initDomain.x1--;
    initDomain.y0++;
    initDomain.y1--;
    initDomain.z0++;
    initDomain.z1--;



    MultiScalarField3D<T> volfracField(nx, ny, nz);
    MultiScalarField3D<T> densityField(nx, ny, nz);
    MultiScalarField3D<T> T_t(volfracField), T_tp1(volfracField), Q(volfracField), v_sed(volfracField);
    MultiScalarField3D<T> D_t(densityField), D_tp1(densityField), Q_d(densityField);
    MultiTensorField3D<T, 3> velocity(volfracField);

    Box3D bottom(0,nx-1,0,ny-1,0,0);
    Box3D top(0,nx-1,0,ny-1,nz-1,nz-1);

    Box3D front(nx-1,nx-1,0,ny-1,1,nz-2);
    Box3D back(0,0,0,ny-1,1,nz-2);

    Box3D left(0,nx-1,0,0,1,nz-2);
    Box3D right(0,nx-1,ny-1,ny-1,1,nz-2);

    nsLattice.toggleInternalStatistics(false);

    ExpSetup(nsLattice, volfracField, densityField, parameters);

    Dynamics<T,NSDESCRIPTOR>* nsbbDynamics = new BounceBack<T,NSDESCRIPTOR>(rho0f);
    defineDynamics(nsLattice, bottom, nsbbDynamics->clone());
    defineDynamics(nsLattice, top, nsbbDynamics->clone());
    defineDynamics(nsLattice, left, nsbbDynamics->clone());
    defineDynamics(nsLattice, right, nsbbDynamics->clone());
    defineDynamics(nsLattice, front, nsbbDynamics->clone());
    defineDynamics(nsLattice, back, nsbbDynamics->clone());
    delete nsbbDynamics; nsbbDynamics = 0;


    setToConstant(volfracField, bottom, 0.);
    setToConstant(volfracField, top, 0.);
    setToConstant(volfracField, left, 0.);
    setToConstant(volfracField, right, 0.);
    setToConstant(volfracField, back, 0.);
    setToConstant(volfracField, front, 0.);

    setToConstant(densityField, bottom, 0.);
    setToConstant(densityField, top, 0.);
    setToConstant(densityField, left, 0.);
    setToConstant(densityField, right, 0.);
    setToConstant(densityField, back, 0.);
    setToConstant(densityField, front, 0.);

      Array<T,NSDESCRIPTOR<T>::d> forceOrientation(T(),T(),(T)1);




    integrateProcessingFunctional (
            new ScalarBuoyanTermProcessor3D<T,NSDESCRIPTOR> (
                tunedg, rho0,
                rhoP, forceOrientation),
            initDomain,
            nsLattice, volfracField, 1 );




    T tIni = global::timer("simTime").stop();
    pcout << "time elapsed for ExpSetup:" << tIni << endl;
    global::timer("simTime").start();

    plint evalTime =5000;
    plint iT = 0;
    plint spongeiT = 1/parameters.getDeltaT();
    plint maxT = 20/parameters.getDeltaT();
    plint saveIter = 0.05/parameters.getDeltaT();;
    util::ValueTracer<T> converge((T)1,(T)100,1.0e-3);

    pcout << "Max Number of iterations: " << maxT << endl;
    pcout << "Number of saving iterations: " << saveIter << endl;


    for (iT = 0; iT <= maxT; ++iT)
    {

        if (iT == (evalTime))
        {
            T tEval = global::timer("simTime").stop();
            T remainTime = (tEval - tIni) / (T)evalTime * (T)maxT/(T)3600;
            global::timer("simTime").start();
            pcout << "Remaining " << (plint)remainTime << " hours, and ";
            pcout << (plint)((T)60*(remainTime - (T)((plint)remainTime))+0.5) << " minutes." << endl;
        }




        if (iT % saveIter == 0)
        {


            //pcout << iT * parameters.getDeltaT() << " : Writing VTK." << endl;
            //writeVTK(nsLattice, volfracField, densityField, parameters, iT);


            pcout << iT * parameters.getDeltaT() << " : Writing gif." << endl;
            writeGif(nsLattice,volfracField,densityField, v_sed, iT);



        }

        if(iT <= spongeiT) omega=0.3;
        else
        omega = nsOmega;

            T spongeZoneFraction = 0.1;
            plint numOutletSpongeCells = util::roundToInt(spongeZoneFraction*nz);
            if (numOutletSpongeCells > 0) {
                // Number of sponge zone lattice nodes at all the outer domain boundaries.
                // So: 0 means the boundary at x = 0
                //     1 means the boundary at x = nx-1
                //     2 means the boundary at y = 0
                //     and so on...


                Array<plint,6> numSpongeCells;
                numSpongeCells[0] = 0;
                numSpongeCells[1] = 0;
                numSpongeCells[2] = 0;
                numSpongeCells[3] = 0;
                numSpongeCells[4] = numOutletSpongeCells;
                numSpongeCells[5] = 0;

                std::vector<MultiBlock3D*> args;
                args.push_back(&nsLattice);

                applyProcessingFunctional (
                    new ViscositySpongeZone3D<T,NSDESCRIPTOR> (
                        nx, ny, nz, omega, numSpongeCells ),
                    nsLattice.getBoundingBox(), args );

            }




        // Lattice Boltzmann iteration step.
        nsLattice.collideAndStream();

        computeVelocity(nsLattice, velocity, initDomain);

        plb::copy(densityField, D_t, initDomain);
        plint numIter =1;
        for (plint iter=0; iter<numIter; ++iter) {
          plb::copy(densityField, D_tp1, initDomain);
          std::vector<MultiBlock3D*> args_d;
    	    args_d.push_back(&D_t);
    	    args_d.push_back(&D_tp1);
    	    args_d.push_back(&densityField);
    	    args_d.push_back(&velocity);
    	    args_d.push_back(&Q_d);
            bool upwind = true;
    	    applyProcessingFunctional(new AdvectionDiffusionFd3D<T>(parameters.getLatticeKappa(), upwind),
                    initDomain, args_d);

          setToConstant(densityField, bottom, 0.);
          setToConstant(densityField, top, 0.);
          setToConstant(densityField, left, 0.);
          setToConstant(densityField, right, 0.);
          setToConstant(densityField, back, 0.);
          setToConstant(densityField, front, 0.);
        }




        plb::copy(volfracField, v_sed, initDomain);
        std::vector<MultiBlock3D*> coupling_args_1;
        coupling_args_1.push_back(&v_sed);
        coupling_args_1.push_back(&densityField);
        applyProcessingFunctional(new DensVsedProcessor3D<T, T>(rhoP, convers), initDomain, coupling_args_1);
        std::vector<MultiBlock3D*> coupling_args_2;
        coupling_args_2.push_back(&v_sed);
        coupling_args_2.push_back(&velocity);
        applyProcessingFunctional(new Velocity_CouplingProcessor3D<T, T, 3>(), initDomain, coupling_args_2);

        plb::copy(volfracField, T_t, initDomain);
        for (plint iter=0; iter<numIter; ++iter) {
          plb::copy(volfracField, T_tp1, initDomain);
          std::vector<MultiBlock3D*> args;
    	    args.push_back(&T_t);
    	    args.push_back(&T_tp1);
    	    args.push_back(&volfracField);
    	    args.push_back(&velocity);
    	    args.push_back(&Q);
            bool upwind = true;
    	    applyProcessingFunctional(new AdvectionDiffusionFd3D<T>(parameters.getLatticeKappa(), upwind),
                    initDomain, args);

          setToConstant(volfracField, bottom, 0.);
          setToConstant(volfracField, top, 0.);
          setToConstant(volfracField, left, 0.);
          setToConstant(volfracField, right, 0.);
          setToConstant(volfracField, back, 0.);
          setToConstant(volfracField, front, 0.);


        }


    }

    writeGif(nsLattice,volfracField,densityField, v_sed, iT);

    T tEnd = global::timer("simTime").stop();

    T totalTime = tEnd-tIni;
    T nx100 = nsLattice.getNx()/(T)100;
    T ny100 = nsLattice.getNy()/(T)100;
    T nz100 = nsLattice.getNz()/(T)100;
    pcout << "N=" << resolution << endl;
    pcout << "number of processors: " << global::mpi().getSize() << endl;
    pcout << "simulation time: " << totalTime << endl;
    pcout << "total time: " << tEnd << endl;
    pcout << "total iterations: " << iT << endl;
    pcout << "Msus: " << nx100*ny100*nz100*(T)iT/totalTime << endl;

    return 0;
}
