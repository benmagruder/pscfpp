#ifndef RPC_SYSTEM_H
#define RPC_SYSTEM_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2015 - 2025, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

// Header file includes 
#include <util/param/ParamComposite.h>    // base class

#include <rpc/solvers/Mixture.h>          // member
#include <rpc/field/Domain.h>             // member
#include <rpc/field/WFieldContainer.h>    // member
#include <rpc/field/CFieldContainer.h>    // member
#include <rpc/field/Mask.h>               // member

#include <prdc/cpu/RField.h>              // member (tmpFieldsRGrid_)
#include <prdc/cpu/RFieldDft.h>           // member (tmpFieldsKGrid_)

#include <pscf/chem/PolymerModel.h>       // member (polymerModel_)

#include <util/misc/FileMaster.h>         // member
#include <util/containers/DArray.h>       // member (tmpFields...)

// Forward references 
namespace Util {
   template <typename T, int N> class FSArray;
}
namespace Pscf {
   class Interaction;
   namespace Prdc {
      template <int D> class UnitCell;
   }
   namespace Rpc {
      template <int D> class Iterator;
      template <int D> class IteratorFactory;
      template <int D> class Sweep;
      template <int D> class SweepFactory;
      template <int D> class Simulator;
      template <int D> class SimulatorFactory;
   }
}

namespace Pscf {
namespace Rpc {

   // Namespaces that are implicitly available, without qualification
   using namespace Util;
   using namespace Pscf;
   using namespace Prdc;
   using namespace Prdc::Cpu;

   /**
   * Main class for SCFT or PS-FTS simulation of one system.
   *
   * A System has (among other member components):
   *
   *    - a Mixture (container for polymer and solvent solvers)
   *    - an Interaction (list of binary chi parameters)
   *    - a Domain (description of unit cell and discretization)
   *    - a container of monomer chemical potential fields
   *    - a container of monomer concentration fields
   *
   * A System may also optionally own Iterator, Sweep, and Simulator
   * (BdSimulator or McSimulator) components. Iterator and Sweep objects 
   * are only used for SCFT calculations.  A Simulator object is only 
   * used for PS-FTS calculations (i.e., field theoretic simulations that 
   * use a partial saddle-point approximation).
   *
   * System is a class template with an integer template parameter 
   * D = 1, 2, or 3 that represents the dimension of space. Many related
   * components are also class templates of the same type. Names such as
   * System, Mixture, Domain, etc. mentioned above are thus names of 
   * class templates, whereas actual class names are of the form
   * Mixture\<D\>, Interaction\<D\>, etc. with D=1, 2, or 3.
   *
   * Usage of a System\<D\> object in the main program looks something 
   * like this:
   * \code
   *    System<D> system;
   *    system.setOptions(argc, argv);
   *    system.readParam();
   *    system.readCommands();
   * \endcode
   * where argc, and argv are parameters containing information about 
   * command line arguments that must be passed from the main program.
   * This is implemented as function template Pscf::Rpc::run in the
   * file src/rpc/pscf_pc.cpp.
   *
   * See also:
   * <ul>
   *  <li> \ref scft_param_pc_page    "Parameter File: SCFT" </li>
   *  <li> \ref psfts_param_page      "Parameter File: PS-FTS" </li>
   *  <li> \ref rpc_System_page       "Parameter File: Full Format" </li>
   *  <li> \ref scft_command_pc_page  "Command File Format" </li>
   * </ul>
   *
   * \ingroup Pscf_Rpc_Module
   */
   template <int D>
   class System : public ParamComposite
   {

   public:

      /// \name Construction and Destruction
      ///@{

      /**
      * Constructor.
      */
      System();

      /**
      * Destructor.
      */
      ~System();

      ///@}
      /// \name Lifetime (Main Actions)
      ///@{

      /**
      * Process command line options.
      *
      * This function takes the same arguments as any C/C++ main program
      * function. The arguments of the main function should d be passed
      * to this function unaltered, to allow this function to process the
      * command line options.
      *
      * \param argc number of command line arguments
      * \param argv array of pointers to command line arguments
      */
      void setOptions(int argc, char **argv);

      /**
      * Read input parameters (with opening and closing lines).
      *
      * \param in input parameter stream
      */
      virtual void readParam(std::istream& in);

      /**
      * Read input parameters from default param file.
      *
      * This function reads the parameter file set by the -p command
      * line option.
      */
      void readParam();

      /**
      * Read body of parameter block (without opening and closing lines).
      *
      * \param in input parameter stream
      */
      virtual void readParameters(std::istream& in);

      /**
      * Read and process commands from an input stream.
      *
      * \param in command script input stream
      */
      void readCommands(std::istream& in);

      /**
      * Read and process commands from the default command file.
      *
      * This function reads the parameter file set by the -c command
      * line option.
      */
      void readCommands();

      ///@}
      /// \name W Field (Chemical Potential Field) Modifiers
      ///@{

      /**
      * Read chemical potential fields in symmetry adapted basis format.
      *
      * This function opens and reads the file named by the "filename"
      * string parameter, which must contain chemical potential fields
      * in symmetry-adapted basis format. The function sets the system
      * w fields equal to those given in this file, by copying the
      * coefficients of the basis expansion and computing values on a
      * real-space grid (r-grid format). 
      *
      * Upon exit, both w().basis() and w().rgrid() are set, w().hasData 
      * and w().isSymmetric() return true, while hasCFields() and 
      * hasFreeEnergy() return false. Unit cell parameters are set to 
      * values read from the field file header.
      * 
      * SCFT calculations that use an iterator that preserves space group
      * symmetry must set an initial field using a function that creates
      * fields that can be represented in symmetry adapted basis form,
      * such as this function, setWBasis, or estimateWFromC.
      *
      * \param filename name of input w-field file in basis format
      */
      void readWBasis(const std::string & filename);

      /**
      * Read chemical potential fields in real space grid (r-grid) format.
      *
      * This function opens and reads the file with the name given by the
      * "filename" string, which must contain chemical potential fields
      * in real space grid (r-grid) format. The function sets values for
      * the system w fields in r-grid format. It does not attempt to set
      * field values in symmetry-adapted basis format, because it cannot
      * assume that the r-grid field exhibits the declared space group
      * symmetry.  On exit, w().rgrid() is reset and w().hasData() 
      * returns true, while w().isSymmetric(), hasCFields() and 
      * hasFreeEnergy() return false. Unit cell parameters are set to 
      * values read from the field file header.
      *
      * Initial chemical potential fields for field theoretic simulations 
      * are normally initialized using a function that sets the fields
      * in r-grid format, such as this function or setWRGrid.
      *
      * \param filename  name of input w-field file in r-grid format
      */
      void readWRGrid(const std::string & filename);

      /**
      * Set chemical potential fields, in symmetry-adapted basis format.
      *
      * This function sets values for w fields in both symmetry adapted
      * and r-grid format by copying coefficient values provided in the 
      * "fields" container that is passed as an argument, and computing
      * values on a real-space grid. Upon return, values of both 
      * w().basis() and w().rgrid() are set, while w().hasData() and 
      * w().isSymmetric() return true, and hasCFields() and hasFreeEnergy()
      * hasFreeEnergy() return false. Unit cell parameters are left 
      * unchanged.
      *
      * \param fields  array of new w (chemical potential) fields
      */
      void setWBasis(DArray< DArray<double> > const & fields);

      /**
      * Set new w fields, in real-space (r-grid) format.
      *
      * This function set values for w fields in r-grid format, but does
      * not set components the symmetry-adapted basis format. Upon return,
      * w.rgrid() is reset and w().hasData() returns true, while
      * w().isSymmetric(), hasCFields() and hasFreeEnergy() all return
      * false.  Unit cell parameters are left unchanged.
      *
      * \param fields  array of new w (chemical potential) fields
      */
      void setWRGrid(DArray< RField<D> > const & fields);

      /**
      * Construct trial w-fields from c-fields in symmetry-adapted form.
      *
      * This function reads concentration fields in symmetrized basis
      * format and constructs an initial guess for corresponding chemical
      * potential fields by setting the Lagrange multiplier field xi to
      * zero. The result is stored in the System w fields container.
      *
      * Upon return, w().hasData() and w().isSymmetric() return true, 
      * while hasCFields() and hasFreeEnergy() return false. Unit
      * cell parameters are set to those read from the c field file
      * header.
      *
      * \param filename  name of input c-field file (basis format)
      */
      void estimateWfromC(const std::string& filename);

      ///@}
      /// \name Unit Cell Modifiers
      ///@{

      /**
      * Set parameters of the associated unit cell.
      *
      * The lattice (i.e., lattice system type) set in the UnitCell<D>
      * unitCell input parameter must agree with any lattice enum value 
      * that was set previously in the parameter file, or an Exception 
      * is thrown.
      *
      * If a space group has been set but a basis has not yet been
      * constructed, then this and the other setUnitCell member functions 
      * wll all construct a symmetry-adapted basis and then allocate
      * memory for fields stored in a symmetry adapted basis format.
      *
      * \param unitCell  new UnitCell<D> (i.e., new parameters)
      */
      void setUnitCell(UnitCell<D> const & unitCell);

      /**
      * Set state of the associated unit cell.
      *
      * The "lattice" (lattice system) enumeration parameter must agree 
      * with any lattice value that was set previously in the parameter 
      * file.
      *
      * See description of setUnitCell(UnitCell<D> const &) for other
      * side effects.
      *
      * \param lattice  lattice system
      * \param parameters  array of new unit cell parameters.
      */
      void setUnitCell(typename UnitCell<D>::LatticeSystem lattice,
                       FSArray<double, 6> const & parameters);

      /**
      * Set parameters of the associated unit cell.
      *
      * The lattice type must have been set before this function is
      * called.  The logical size of the FSArray<double, 6> "parameters" 
      * array must match the expected number of parameters for the 
      * current lattice type.
      *
      * See description of setUnitCell(UnitCell<D> const &) for other
      * side effects.
      *
      * \param parameters  array of new unit cell parameters.
      */
      void setUnitCell(FSArray<double, 6> const & parameters);

      ///@}
      /// \name Primary Field Theory Computations
      ///@{

      /**
      * Solve the modified diffusion equation once, without iteration.
      *
      * This function calls the Mixture::compute() function to solve
      * the statistical mechanics problem for a non-interacting system
      * subjected to the currrent chemical potential fields. This
      * requires solution of the modified diffusion equation for all
      * polymers, computation of Boltzmann weights for all solvents,
      * computation of molecular partition functions for all species,
      * computation of concentration fields for blocks and solvents,
      * and computation of overall concentrations for all monomer types.
      * This function does not compute the canonical (Helmholtz) free
      * energy or grand-canonical free energy (i.e., pressure). Upon
      * return, the flag hasCFields is set true.
      *
      * If argument needStress == true, then this function also calls
      * Mixture<D>::computeStress() to compute the stress.
      *
      * \pre The w().hasData() flag must be true on entry, to confirm
      * that chemical potential fields have been set.
      *
      * \param needStress true if stress is needed, false otherwise
      */
      void compute(bool needStress = false);

      /**
      * Iteratively solve a SCFT problem.
      *
      * This function calls the iterator to solve the SCFT problem for
      * the current mixture and system parameters, using the current
      * chemical potential fields and current unit cell parameter values 
      * as initial guesses.  On exit, hasCFields is set true whether or 
      * not convergence is obtained to within the desired tolerance.  
      * The Helmholtz free energy and pressure are computed only if
      * convergence is obtained.
      *
      * \pre Function hasIterator() must return true.
      * \pre Function w().hasData() flag must return true.
      * \pre Function w().isSymmetric() flag must return true if the chosen
      * iterator uses a symmetry adapted basis, and so requires this.
      *
      * \return returns 0 for successful convergence, 1 for failure.
      *
      * \param isContinuation true if continuation within a sweep.
      */
      int iterate(bool isContinuation = false);

      /**
      * Sweep in parameter space, solving an SCF problem at each point.
      *
      * This function uses a Sweep object that was initialized in the
      * parameter file to solve the SCF problem at a sequence of points
      * along a contour in parameter space. The nature of this sequence
      * of points is determined by implementation of a subclass of Sweep
      * and the parameters passed to the sweep object in the parameter
      * file.  The Iterator that is initialized in the parameter file
      * is called at each state point.
      *
      * \pre Function hasSweep() must return true.
      * \pre All preconditions of the iterate() function must be satisfied
      */
      void sweep();

      /**
      * Perform a field theoretic simulation (PS-FTS).
      *
      * Perform a field theoretic simulation using the partial saddle-point
      * approximation (PS-FTS). The type of calculation (BD or MC) is
      * determined by the type of Simulator (BdSimulator or McSimulator)
      * that created in the parameter file. The number of BD steps or 
      * attempted MC moves to be performed is given by the parameter
      * "nStep" .
      *
      * \pre Function hasSimulator() must return true. 
      * \pre Function w().hasData() must return true. 
      *
      * \param nStep  number of simulation steps
      */
      void simulate(int nStep);

      ///@}
      /// \name Thermodynamic Properties
      ///@{

      /**
      * Compute free energy density and pressure for current fields.
      *
      * This function should be called after a successful call of
      * System::iterate() or Iterator::solve(). Resulting values of
      * Helmholtz free energy per monomer, free energy components fIdeal,
      * fInter, and fExt, and the non-dimensional pressure (-1 times the
      * grand canonical free energy per monomer) are stored in private
      * member variables and can be retrieved by the member functions
      * fHelmholtz(), fIdeal(), fInter(), fExt(), and  pressure().
      */
      void computeFreeEnergy();

      /**
      * Get precomputed Helmoltz free energy per monomer / kT.
      *
      * The value retrieved by this function is computed by the
      * computeFreeEnergy() function.
      */
      double fHelmholtz() const;

      /**
      * Get the ideal gas contribution to fHelmholtz(). 
      *
      * This function retrieves a value computed by computeFreeEnergy().
      */
      double fIdeal() const;

      /**
      * Get the interaction contribution to fHelmholtz(). 
      *
      * This function retrieves a value computed by computeFreeEnergy().
      */
      double fInter() const;

      /**
      * Get the external field contribution to fHelmholtz(). 
      *
      * This function retrieves a value computed by computeFreeEnergy().
      */
      double fExt() const;

      /**
      * Get precomputed pressure x monomer volume / kT.
      *
      * The value retrieved by this function is computed by the
      * computeFreeEnergy() function. The value is -1 times the grand
      * canonical free energy per monomer divided by kT. 
      */
      double pressure() const;

      ///@}
      /// \name Thermodynamic Data Output
      ///@{

      /**
      * Write parameter file to an ostream, omitting any sweep block.
      *
      * This function writes the Mixture, Interaction, and Domain blocks 
      * of a parameter file, as well as any Iterator block, but omits any
      * Sweep or Simulator blocks. The intent is to produce an output 
      * during an SCFT sweep that only refers to parameters relevant to a 
      * single state point, in form that could be used a a parameter file 
      * for a a single SCFT calculation.
      *
      * \param out output stream
      */
      void writeParamNoSweep(std::ostream& out) const;

      /**
      * Write thermodynamic properties to a file.
      *
      * This function outputs Helmholtz free energy per monomer, pressure
      * (in units of kT per monomer volume), the volume fraction and
      * chemical potential of each species, and all unit cell parameters.
      *
      * Calling writeParamNoSweep and writeThermo in succession with the 
      * same output stream will produce a single file containing both 
      * input parameters and resulting thermodynamic properties.
      *
      * \param out output stream
      */
      void writeThermo(std::ostream& out);
      
      /**
      * Write stress properties to a file.
      *
      * This function outputs derivatives of SCFT free energy w/ respect to 
      * each unit cell parameters.
      *
      * The log output created by an Iterator upon convergence should call
      * writeStress after writeThermo if and only if the iterator is not 
      * flexible. 
      *
      * \param out output stream
      */
      void writeStress(std::ostream& out);

      ///@}
      /// \name Field Output
      ///@{

      /**
      * Write chemical potential fields in symmetrized basis format.
      *
      * \param filename name of output file
      */
      void writeWBasis(const std::string & filename) const;

      /**
      * Write chemical potential fields in real space grid (r-grid) format.
      *
      * \param filename name of output file
      */
      void writeWRGrid(const std::string & filename) const;

      /**
      * Write concentration fields in symmetrized basis format.
      *
      * \param filename name of output file
      */
      void writeCBasis(const std::string & filename) const;

      /**
      * Write concentration fields in real space grid (r-grid) format.
      *
      * \param filename name of output file
      */
      void writeCRGrid(const std::string & filename) const;

      /**
      * Write c-fields for all blocks and solvents in r-grid format.
      *
      * Writes concentrations for all blocks of all polymers and all
      * solvent species in r-grid format. Columns associated with blocks
      * appear ordered by polymer id and block id, with blocks of the
      * same polymer listed sequentially, followed by columns associated
      * with solvent species ordered by solvent id.
      *
      * \param filename name of output file
      */
      void writeBlockCRGrid(const std::string & filename) const;

      ///@}
      /// \name Propagator Output
      ///@{

      /**
      * Write slice of a propagator at fixed s in r-grid format.
      *
      * \param filename  name of output file
      * \param polymerId  integer id of the polymer
      * \param blockId  integer id of the block within the polymer
      * \param directionId  integer id of the direction (0 or 1)
      * \param segmentId  integer integration step index
      */
      void writeQSlice(std::string const & filename,
                       int polymerId, int blockId,
                       int directionId, int segmentId)  const;

      /**
      * Write the final slice of a propagator in r-grid format.
      *
      * \param filename  name of output file
      * \param polymerId  integer id of the polymer
      * \param blockId  integer id of the block within the polymer
      * \param directionId  integer id of the direction (0 or 1)
      */
      void writeQTail(std::string const & filename, int polymerId,
                      int blockId, int directionId)  const;

      /**
      * Write one propagator for one block, in r-grid format.
      *
      * \param filename  name of output file
      * \param polymerId  integer id of the polymer
      * \param blockId  integer id of the block within the polymer
      * \param directionId  integer id of the direction (0 or 1)
      */
      void writeQ(std::string const & filename, int polymerId,
                  int blockId, int directionId)  const;

      /**
      * Write all propagators of all blocks, each to a separate file.
      *
      * Write all propagators for both directions for all blocks
      * of all polymers, with each propagator in a separate file.
      * The function writeQ is called internally for each propagator,
      * and is passed an automatically generated file name.  The file
      * name for each propagator is given by a string of the form
      * (basename)_(ip)_(ib)_(id), where (basename) denotes the value
      * of the std::string function parameter basename, and where
      * (ip), (ib), and (id) denote the string representations of
      * a polymer indiex ip, a block index ib, and direction index id,
      * with id = 0 or 1. For example, if basename == "out/q", then
      * the file name of the propagator for direction 1 of block 2
      * of polymer 0 would be "out/q_0_2_1".
      *
      * \param basename  common prefix for output file names
      */
      void writeQAll(std::string const & basename);

      ///@}
      /// \name Field File Manipulation and Analysis
      ///@{

      /**
      * Convert a field from symmetrized basis format to r-grid format.
      *
      * This function reads a field file in basis format, converts the
      * fields to r-grid format, and writes the fields in r-grid format
      * to a different file.
      *
      * This and other field conversion functions do not change the w
      * or c fields stored by this System - all required calculations
      * are performed using mutable temporary memory.
      *
      * \param inFileName name of input file (basis format)
      * \param outFileName name of output file (r-grid format)
      */
      void basisToRGrid(const std::string & inFileName,
                        const std::string & outFileName);

      /**
      * Convert a field from real-space grid to symmetrized basis format.
      *
      * This function checks if the input fields have the declared space
      * group symmetry, and prints a warning if it detects deviations
      * that exceed some small threshhold, but proceeds to attempt the
      * conversion even if such an error is detected. Converting a field
      * that does not have the declared space group symmetry to basis
      * format is a destructive operation that modifies the field in
      * unpredictable ways.
      *
      * \param inFileName name of input file (r-grid format)
      * \param outFileName name of output file (basis format)
      */
      void rGridToBasis(const std::string & inFileName,
                        const std::string & outFileName);

      /**
      * Convert fields from Fourier (k-grid) to real-space (r-grid) format.
      *
      * \param inFileName name of input file (k-grid format)
      * \param outFileName name of output file (r-grid format)
      */
      void kGridToRGrid(const std::string& inFileName,
                        const std::string& outFileName);

      /**
      * Convert fields from real-space (r-grid) to Fourier (k-grid) format.
      *
      * \param inFileName name of input file (r-grid format)
      * \param outFileName name of output file (k-grid format)
      */
      void rGridToKGrid(const std::string & inFileName,
                        const std::string & outFileName);

      /**
      * Convert fields from Fourier (k-grid) to symmetrized basis format.
      *
      * This function checks if the input fields have the declared space
      * group symmetry, and prints a warning if it detects deviations
      * that exceed some small threshhold, but proceeds to attempt the
      * conversion even if such an error is detected. Converting a field
      * that does not have the declared space group symmetry to basis
      * format is a destructive operation that modifies the field in
      * unpredictable ways.
      *
      * \param inFileName name of input file (k-grid format)
      * \param outFileName name of output file (basis format)
      */
      void kGridToBasis(const std::string& inFileName,
                        const std::string& outFileName);

      /**
      * Convert fields from symmetrized basis to Fourier (k-grid) format.
      *
      * \param inFileName name of input file (basis format)
      * \param outFileName name of output file (k-grid format)
      */
      void basisToKGrid(const std::string & inFileName,
                        const std::string & outFileName);

      /**
      * Compare arrays of fields in basis format, output a report.
      *
      * Outputs maximum and root-mean-squared differences to the
      * standard Log file. 
      *
      * \param field1  first array of fields (basis format)
      * \param field2  second array of fields (basis format)
      */
      void compare(const DArray< DArray<double> > field1,
                   const DArray< DArray<double> > field2);

      /**
      * Compare two fields in r-grid format, output a report.
      *
      * Outputs maximum and root-mean-squared differences to the
      * standard Log file. 
      *
      * \param field1  first array of fields (r-grid format)
      * \param field2  second array of fields (r-grid format)
      */
      void compare(const DArray< RField<D> > field1,
                   const DArray< RField<D> > field2);

      /**
      * Check if r-grid fields have the declared space group symmetry.
      *
      * \param inFileName name of input file
      * \param epsilon error threshold used when testing for symmetry
      * \return true if fields all have symmetry, false otherwise
      */
      bool checkRGridFieldSymmetry(const std::string & inFileName,
                                   double epsilon = 1.0E-8);

      /**
      * Multiply all components of an array of basis fields by a scalar.
      *
      * \param inFileName  name of input field file 
      * \param outFileName  name of file for rescaled output fields
      * \param factor  factor by which to multiply all field elements
      */
      void scaleFieldsBasis(const std::string & inFileName,
                            const std::string & outFileName,
                            double factor);

      /**
      * Multiply all elements of an array of r-grid fields by a scalar.
      *
      * \param inFileName  name of input field file 
      * \param outFileName  name of file for rescaled output fields
      * \param factor  factor by which to multiply all field elements
      */
      void scaleFieldsRGrid(const std::string & inFileName,
                            const std::string & outFileName,
                            double factor) const;

      /**
      * Expand the number of spatial dimensions of an r-grid field.
      *
      * This function reads a D-dimensional field and outputs a field
      * in a format appropriate for d-dimensional space, for d > D, by
      * assuming that all field values are independent of coordinates
      * associated with the added dimensions. It can thus create a file
      * representing a field with lamellar (D=1) or hexagonal (D=2)
      * symmetry on a 3D (d=3) grid.
      *
      * Element i of array newGridDimensions contains the number of
      * grid points in added dimension D + i. This array must have a
      * capacity d - D.
      *
      * \param inFileName filename name of input field file
      * \param outFileName filename name of output field file
      * \param d  intended dimensions (d > D)
      * \param newGridDimensions number of grid points in added dimensions
      */
      void expandRGridDimension(const std::string & inFileName,
                                const std::string & outFileName,
                                int d,
                                DArray<int> newGridDimensions);
      /**
      * Replicate the crystal unit cell to create a larger cell. 
      *
      * This function reads a D-dimensional field and replicates the
      * unit cell a specified number of times in each D direction
      *
      * Element i of array replicas contains the number of replication
      * times in direction i.
      *
      * \param inFileName filename name of input field file
      * \param outFileName filename name of output field file
      * \param replicas  the number of replicas in each D direction
      */
      void replicateUnitCell(const std::string & inFileName,
                             const std::string & outFileName,
                             IntVec<D> const & replicas);

      ///@}
      /// \name Timers
      ///@{

      /**
      * Write timer information to an output stream.
      *
      * \param out output stream
      */
      void writeTimers(std::ostream& out);

      /**
      * Clear timers
      */
      void clearTimers();

      ///@}
      /// \name Field Accessors
      ///@{

      /**
      * Get container of chemical potential fields (w fields).
      */
      WFieldContainer<D> const & w() const;

      /**
      * Get container of monomer concentration fields (c fields).
      */
      CFieldContainer<D> const & c() const;

      /**
      * Get container of external potential fields (reference).
      */
      WFieldContainer<D>& h();

      /**
      * Get container of external potential fields (const reference).
      */
      WFieldContainer<D> const & h() const;

      /**
      * Get the mask (field to which total density is constrained).
      */
      Mask<D>& mask();

      /**
      * Get the mask by const reference.
      */
      Mask<D> const & mask() const;

      ///@}
      /// \name Member Object Accessors
      ///@{

      /**
      * Get the Mixture by non-const reference.
      */
      Mixture<D> & mixture();

      /**
      * Get the Mixture by const reference.
      */
      Mixture<D> const & mixture() const;

      /**
      * Get Interaction (excess free energy model) by reference.
      */
      Interaction& interaction();

      /**
      * Get Interaction (excess free energy model) by const reference.
      */
      Interaction const & interaction() const;

      /**
      * Get Domain by const reference.
      */
      Domain<D> const & domain() const;

      /**
      * Get the Iterator by non-const reference.
      */
      Iterator<D>& iterator();

      /**
      * Get the Iterator by const reference.
      */
      Iterator<D> const & iterator() const;

      /**
      * Get the Simulator for field theoretic simulation.
      */
      Simulator<D>& simulator();

      /**
      * Get the FileMaster by non-const references.
      *
      * Access by non-const reference is used in some unit tests.
      */
      FileMaster& fileMaster();

      /**
      * Get the FileMaster by const reference.
      */
      FileMaster const & fileMaster() const;

      ///@}
      /// \name Boolean Queries
      ///@{

      /**
      * Does this system have an Iterator object?
      */
      bool hasIterator() const;

      /**
      * Does this system have a Sweep object?
      */
      bool hasSweep() const;

      /**
      * Does this system have an initialized Simulator?
      */
      bool hasSimulator() const;

      /**
      * Does this system have external potential fields?
      */
      bool hasExternalFields() const;

      /**
      * Does this system have a mask (inhomogeneous density constraint)?
      */
      bool hasMask() const;

      /**
      * Have c fields been computed from the current w fields?
      */
      bool hasCFields() const;

      /**
      * Has the free energy been computed from the current w fields?
      */
      bool hasFreeEnergy() const;

      ///@}

   private:

      // Component objects

      /**
      * Mixture object (solves MDE for all species).
      */
      Mixture<D> mixture_;

      /**
      * Domain object (unit cell, space group, mesh, and basis).
      */
      Domain<D> domain_;

      /**
      * Filemaster (holds paths to associated I/O files).
      */
      FileMaster fileMaster_;

      // Pointers to objects owned by the System

      /**
      * Pointer to Interaction (free energy model).
      */
      Interaction* interactionPtr_;

      /**
      * Pointer to an iterator.
      */
      Iterator<D>* iteratorPtr_;

      /**
      * Pointer to iterator factory object
      */
      IteratorFactory<D>* iteratorFactoryPtr_;

      /**
      * Pointer to a Sweep object
      */
      Sweep<D>* sweepPtr_;

      /**
      * Pointer to SweepFactory object
      */
      SweepFactory<D>* sweepFactoryPtr_;

      /**
      * Simulator - coordinator for field theoretic simulation.
      */
      Simulator<D>* simulatorPtr_;

      /**
      * Pointer to Simulator factory object
      */
      SimulatorFactory<D>* simulatorFactoryPtr_;

      // Field container objects

      /**
      * Chemical potential fields.
      */
      WFieldContainer<D> w_;

      /**
      * Monomer concentration / volume fraction fields.
      */
      CFieldContainer<D> c_;

      /**
      * External potential fields.
      */
      WFieldContainer<D> h_;

      /**
      * Field to which the total density is constrained.
      */
      Mask<D> mask_;

      /**
      * Work array of field coefficients for all monomer types.
      *
      * Indexed by monomer typeId, size = nMonomer.
      */
      mutable DArray< DArray<double> > tmpFieldsBasis_;

      /**
      * Work array of fields on real space grid.
      *
      * Indexed by monomer typeId, size = nMonomer.
      */
      mutable DArray< RField<D> > tmpFieldsRGrid_;

      /**
      * Work array of fields on Fourier grid (k-grid).
      *
      * Indexed by monomer typeId, size = nMonomer.
      */
      mutable DArray< RFieldDft<D> > tmpFieldsKGrid_;

      // Thermodynamic properties

      /**
      * Helmholtz free energy per monomer / kT.
      */
      double fHelmholtz_;

      /**
      * Ideal gas contribution to fHelmholtz_.
      *
      * This encompasses the internal energy and entropy of
      * non-interacting free chains in their corresponding
      * potential fields defined by w_.
      */
      double fIdeal_;

      /**
      * Multi-chain interaction contribution to fHelmholtz_.
      */
      double fInter_;

      /**
      * External field contribution to fHelmholtz_ (if any).
      */
      double fExt_;

      /**
      * Pressure times monomer volume / kT.
      *
      * This is -1 times the grand-canonical free energy per monomer,
      * divided by kT.
      */
      double pressure_;

      // Boolean enumerations variables

      /**
      * Has memory been allocated for fields in grid format?
      */
      bool isAllocatedGrid_;

      /**
      * Has memory been allocated for fields in symmetrized basis format?
      */
      bool isAllocatedBasis_;

      /**
      * Has the mixture been initialized?
      */
      bool hasMixture_;

      /**
      * Have c fields been computed for the current w fields?
      *
      * Set this true when c fields are computed by solving the MDEs for
      * for all blocks. Set this false whenever the system w fields in
      * w_ or the unit cell parameters in the Domain are modified. When 
      * hasCFields_ is true, both the c fields for individual blocks and 
      * solvent species in the Mixture and the total fields for different 
      * monomer types in the System::c_ container are those computed from 
      * the current w fields in System::w_ container, using the current
      * unit cell parameters.
      */
      bool hasCFields_;

      /**
      * Has fHelmholtz been computed for the current w and c fields?
      */
      bool hasFreeEnergy_;

      /**
      * Value for polymer model (thread or bead) read from file.
      */
      PolymerModel::Type polymerModel_;

      // Private member functions

      /**
      * Allocate memory for fields in grid formats (private)
      */
      void allocateFieldsGrid();

      /**
      * Allocate memory for fields in basis format (private)
      */
      void allocateFieldsBasis();

      /**
      * Read a field file header, make the basis if not done previously.
      *
      * Used to peek at a file header to get initial unit cell parameters,
      * use this to initialize basis if not done previously.
      *
      * \param filename name of field file
      */
      void readFieldHeader(std::string filename);

      /**
      * Read a string and echo to log file.
      *
      * Used to read filenames in readCommands.
      *
      * \param in  input stream (i.e., input file)
      * \param string  string to read and echo
      */
      void readEcho(std::istream& in, std::string& string) const;

      /**
      * Read a floating point number and echo to log file.
      *
      * Used to read filenames in readCommands.
      *
      * \param in  input stream (i.e., input file)
      * \param value  number to read and echo
      */
      void readEcho(std::istream& in, double& value) const;

   };

   // Inline member functions

   // Get the Mixture object by non-const reference.
   template <int D>
   inline Mixture<D>& System<D>::mixture()
   {  return mixture_; }

   // Get the Mixture object by const reference.
   template <int D>
   inline Mixture<D> const & System<D>::mixture() const
   {  return mixture_; }

   // Get the Domain by const reference.
   template <int D>
   inline Domain<D> const & System<D>::domain() const
   {  return domain_; }

   // Get the Simulator.
   template <int D>
   inline Simulator<D>& System<D>::simulator()
   {  return *simulatorPtr_; }

   // Get the FileMaster by non-const reference.
   template <int D>
   inline FileMaster& System<D>::fileMaster()
   {  return fileMaster_; }

   // Get the FileMaster by const reference.
   template <int D>
   inline FileMaster const & System<D>::fileMaster() const
   {  return fileMaster_; }

   // Get the Interaction (excess free energy model).
   template <int D>
   inline Interaction& System<D>::interaction()
   {
      UTIL_ASSERT(interactionPtr_);
      return *interactionPtr_;
   }

   // Get the Interaction by const reference.
   template <int D>
   inline Interaction const & System<D>::interaction() const
   {
      UTIL_ASSERT(interactionPtr_);
      return *interactionPtr_;
   }

   // Get the Iterator.
   template <int D>
   inline Iterator<D>& System<D>::iterator()
   {
      UTIL_ASSERT(iteratorPtr_);
      return *iteratorPtr_;
   }

   // Get the Iterator by const reference.
   template <int D>
   inline Iterator<D> const & System<D>::iterator() const
   {
      UTIL_ASSERT(iteratorPtr_);
      return *iteratorPtr_;
   }

   // Get container of chemical potential fields (const reference)
   template <int D>
   inline
   WFieldContainer<D> const & System<D>::w() const
   {  return w_; }

   // Get container of monomer concentration fields (const reference)
   template <int D>
   inline
   CFieldContainer<D> const & System<D>::c() const
   {  return c_; }

   // Get container of external potential fields (reference)
   template <int D>
   inline WFieldContainer<D>& System<D>::h()
   {  return h_; }

   // Get container of external potential fields (const reference)
   template <int D>
   inline WFieldContainer<D> const & System<D>::h() const
   {  return h_; }

   // Get mask field (reference)
   template <int D>
   inline Mask<D>& System<D>::mask()
   {  return mask_; }

   // Get mask field (const reference)
   template <int D>
   inline Mask<D> const & System<D>::mask() const
   {  return mask_; }

   // Does the system have an Iterator object?
   template <int D>
   inline bool System<D>::hasIterator() const
   {  return (iteratorPtr_); }

   // Does the system have a Sweep object?
   template <int D>
   inline bool System<D>::hasSweep() const
   {  return (sweepPtr_); }

   // Does this system have external potential fields?
   template <int D>
   inline bool System<D>::hasExternalFields() const
   {  return h_.hasData(); }

   // Does this system have a mask?
   template <int D>
   inline bool System<D>::hasMask() const
   {  return mask_.hasData(); }

   // Does the system have an initialized Simulator ?
   template <int D>
   inline bool System<D>::hasSimulator() const
   {  return (simulatorPtr_); }

   // Have the c fields been computed for the current w fields?
   template <int D>
   inline bool System<D>::hasCFields() const
   {  return hasCFields_; }

   // Get the precomputed Helmholtz free energy per monomer / kT.
   template <int D>
   inline double System<D>::fHelmholtz() const
   {
      UTIL_CHECK(hasFreeEnergy_);
      return fHelmholtz_;
   }

   // Get the ideal gas contribution to fHelmholtz per monomer / kT.
   template <int D>
   inline double System<D>::fIdeal() const
   {
      UTIL_CHECK(hasFreeEnergy_);
      return fIdeal_;
   }

   // Get the interaction contribution to fHelmholtz per monomer / kT.
   template <int D>
   inline double System<D>::fInter() const
   {
      UTIL_CHECK(hasFreeEnergy_);
      return fInter_;
   }

   // Get the external field contribution to fHelmholtz per monomer / kT.
   template <int D>
   inline double System<D>::fExt() const
   {
      UTIL_CHECK(hasFreeEnergy_);
      return fExt_;
   }

   // Get the non-dimensional pressure (units of kT / monomer volume).
   template <int D>
   inline double System<D>::pressure() const
   {
      UTIL_CHECK(hasFreeEnergy_);
      return pressure_;
   }

   // Has the free energy been computed for the current w fields?
   template <int D>
   inline bool System<D>::hasFreeEnergy() const
   {  return hasFreeEnergy_; }

   #ifndef RPC_SYSTEM_TPP
   // Suppress implicit instantiation
   extern template class System<1>;
   extern template class System<2>;
   extern template class System<3>;
   #endif

} // namespace Rpc
} // namespace Pscf
#endif
