#ifndef RPC_FIELD_IO_TPP
#define RPC_FIELD_IO_TPP

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "FieldIo.h"

#include <prdc/cpu/FFT.h>
#include <prdc/cpu/RField.h>
#include <prdc/cpu/RFieldDft.h>

#include <prdc/field/fieldIoUtil.h>
#include <prdc/crystal/fieldHeader.h>

#include <prdc/crystal/shiftToMinimum.h>
#include <prdc/crystal/UnitCell.h>
#include <prdc/crystal/Basis.h>
#include <prdc/crystal/SpaceGroup.h>
#include <prdc/crystal/UnitCell.h>

#include <pscf/mesh/Mesh.h>
#include <pscf/mesh/MeshIterator.h>
#include <pscf/mesh/MeshIteratorFortran.h>
#include <pscf/math/IntVec.h>

#include <util/misc/FileMaster.h>
#include <util/misc/Log.h>
#include <util/format/Str.h>
#include <util/format/Int.h>
#include <util/format/Dbl.h>

#include <iomanip>
#include <string>

namespace Pscf {
namespace Rpc {

   using namespace Util;
   using namespace Pscf::Prdc;
   using namespace Pscf::Prdc::Cpu;

   /*
   * Constructor.
   */
   template <int D>
   FieldIo<D>::FieldIo()
    : meshPtr_(0),
      fftPtr_(0),
      hasGroupPtr_(0),
      groupNamePtr_(0),
      groupPtr_(0),
      basisPtr_(0),
      fileMasterPtr_()
   {}

   /*
   * Destructor.
   */
   template <int D>
   FieldIo<D>::~FieldIo()
   {}

   /*
   * Create associations with other members of parent Domain.
   */
   template <int D>
   void
   FieldIo<D>::associate(
                    Mesh<D> const & mesh,
                    FFT<D> const & fft,
                    typename UnitCell<D>::LatticeSystem const & lattice,
                    bool const & hasGroup,
                    std::string const & groupName,
                    SpaceGroup<D> const & group,
                    Basis<D> & basis)
   {
      meshPtr_ = &mesh;
      fftPtr_ = &fft;
      latticePtr_ = &lattice;
      hasGroupPtr_ = &hasGroup;
      groupNamePtr_ = &groupName;
      groupPtr_ = &group;
      basisPtr_ = &basis;
   }

   /*
   * Create an association with a FileMaster.
   */
   template <int D>
   void FieldIo<D>::setFileMaster(FileMaster const & fileMaster)
   {  fileMasterPtr_ = &fileMaster; }

   /*
   * Read a set of fields in basis format.
   */
   template <int D>
   void FieldIo<D>::readFieldsBasis(std::istream& in,
                                    DArray< DArray<double> >& fields,
                                    UnitCell<D>& unitCell) const
   {
      // Preconditions
      UTIL_CHECK(hasGroup());

      // Read header of field file (also checks group name)
      int nMonomer;
      bool isSymmetric;
      FieldIo<D>::readFieldHeader(in, nMonomer, unitCell, isSymmetric);
      UTIL_CHECK(isSymmetric);
      UTIL_CHECK(basis().isInitialized());

      // Read the number of stars into nStarIn
      std::string label;
      in >> label;
      if (label != "N_star" && label != "N_basis") {
         std::string msg =  "\n";
         msg += "Error reading field file:\n";
         msg += "Expected N_basis or N_star, but found [";
         msg += label;
         msg += "]";
         UTIL_THROW(msg.c_str());
      }
      int nStarIn;
      in >> nStarIn;
      UTIL_CHECK(nStarIn > 0);

      // Check allocation of fields container
      if (fields.isAllocated()) {
         int nMonomerFields = fields.capacity();
         UTIL_CHECK(nMonomerFields > 0);
         UTIL_CHECK(nMonomerFields == nMonomer);
         int fieldCapacity = fields[0].capacity();
         for (int i = 0; i < nMonomer; ++i) {
            UTIL_CHECK( fields[i].capacity() == fieldCapacity );
         }
      } else {
         fields.allocate(nMonomer);
         for (int i = 0; i < nMonomer; ++i) {
            fields[i].allocate(nStarIn);
         }
      }

      // Read data in symmetry-adapted basis format
      readBasisData(in, fields, unitCell, mesh(), basis(), nStarIn);
   }

   /*
   * Open/close a file and read a set of fields in basis format.
   */
   template <int D>
   void FieldIo<D>::readFieldsBasis(std::string filename,
                                    DArray<DArray<double> >& fields,
                                    UnitCell<D>& unitCell) const
   {

      std::ifstream file;
      fileMaster().openInputFile(filename, file);
      readFieldsBasis(file, fields, unitCell);
      file.close();
   }

   /*
   * Read a single fields in basis format into input stream
   */
   template <int D>
   void FieldIo<D>::readFieldBasis(std::istream& in,
                                   DArray<double>& field,
                                   UnitCell<D>& unitCell) const
   {
      // Local array, of data type required by readFieldsBasis
      DArray<DArray<double> > fields;

      // If single field is allocated, allocate local array fields
      if (field.isAllocated()) {
         fields.allocate(1);
         fields[0].allocate(field.capacity());
      }
      // Otherwise, pass unallocated fields array to readFieldsBasis

      // Read file containing a single field, allocate fields if needed.
      readFieldsBasis(in, fields, unitCell);

      // Check that it only read 1 field
      UTIL_CHECK(fields.capacity() == 1);

      // Copy data from local array fields to function parameter field
      field = fields[0];
   }

   /*
   * Open-close a file and read single fields in basis format.
   */
   template <int D>
   void FieldIo<D>::readFieldBasis(std::string filename,
                                   DArray<double>& field,
                                   UnitCell<D>& unitCell)
   const
   {
      std::ifstream file;
      fileMaster().openInputFile(filename, file);
      readFieldBasis(file, field, unitCell);
      file.close();
   }

   /*
   * Write an array of fields in basis format to an output stream.
   */
   template <int D>
   void
   FieldIo<D>::writeFieldsBasis(std::ostream &out,
                                DArray<DArray<double> > const & fields,
                                UnitCell<D> const & unitCell) const
   {
      int nMonomer = fields.capacity();
      int fieldCapacity = fields[0].capacity();
      UTIL_CHECK(nMonomer > 0);
      UTIL_CHECK(fieldCapacity > 0);
      UTIL_CHECK(basis().isInitialized());

      // Write header (common portion)
      bool isSymmetric = true;
      writeFieldHeader(out, nMonomer, unitCell, isSymmetric);

      // Write N_Basis
      int nBasis = basis().nBasis();
      if (nBasis > fieldCapacity) {
         nBasis = fieldCapacity;
      }
      out << "N_basis      " << std::endl
          << "             " << nBasis << std::endl;

      writeBasisData(out, fields, basis());

   }

   /*
   * Open-close a file, and write an array of fields in basis format.
   */
   template <int D>
   void
   FieldIo<D>::writeFieldsBasis(std::string filename,
                                DArray<DArray<double> > const & fields,
                                UnitCell<D> const & unitCell) const
   {
       std::ofstream file;
       fileMaster().openOutputFile(filename, file);
       writeFieldsBasis(file, fields, unitCell);
       file.close();
   }

   /*
   * Write a single field in basis format to an output stream.
   */
   template <int D>
   void FieldIo<D>::writeFieldBasis(std::ostream& out,
                                    DArray<double> const & field,
                                    UnitCell<D> const & unitCell)
   const
   {
      // Create local array of type required by writeFieldsBasis
      DArray<DArray<double> > fields;
      fields.allocate(1);
      fields[0].allocate(field.capacity());

      // Copy data from input parameter to local array 
      fields[0] = field;

      writeFieldsBasis(out, fields, unitCell);
   }

   /*
   * Write a single field in basis format, open and close the file.
   */
   template <int D>
   void
   FieldIo<D>::writeFieldBasis(std::string filename,
                               DArray<double> const & field,
                               UnitCell<D> const & unitCell) const
   {
      // Create local array of type required by writeFieldsBasis
      DArray<DArray<double> > fields;
      fields.allocate(1);
      fields[0].allocate(field.capacity());

      // Copy data from input parameter to local array
      fields[0] = field;

      writeFieldsBasis(filename, fields, unitCell);
   }

   // R-Grid Field Format IO

   /*
   * Read an array of fields in r-grid format from an input stream.
   */
   template <int D>
   void FieldIo<D>::readFieldsRGrid(std::istream &in,
                                    DArray<RField<D> >& fields,
                                    UnitCell<D>& unitCell) const
   {
      int nMonomer;
      bool isSymmetric;
      FieldIo<D>::readFieldHeader(in, nMonomer, unitCell, isSymmetric);
      readMeshDimensions(in);
      checkAllocateFields(fields, mesh().dimensions(), nMonomer);
      readRGridData(in, fields, mesh().dimensions(), nMonomer);

   }

   template <int D>
   void FieldIo<D>::readFieldsRGrid(std::string filename,
                                    DArray< RField<D> >& fields,
                                    UnitCell<D>& unitCell)
   const
   {
      std::ifstream file;
      fileMaster().openInputFile(filename, file);
      readFieldsRGrid(file, fields, unitCell);
      file.close();
   }

   template <int D>
   void FieldIo<D>::readFieldRGrid(std::istream &in,
                                   RField<D> & field,
                                   UnitCell<D>& unitCell)
   const
   {

      // Read header
      int nMonomer;
      bool isSymmetric;
      FieldIo<D>::readFieldHeader(in, nMonomer, unitCell, isSymmetric);
      UTIL_CHECK(nMonomer == 1);
      readMeshDimensions(in);

      // Read data
      checkAllocateField(field, mesh().dimensions());
      readRGridData(in, field, mesh().dimensions());
   }

   template <int D>
   void FieldIo<D>::readFieldRGrid(std::string filename,
                                    RField<D> & field,
                                    UnitCell<D>& unitCell)
   const
   {
      std::ifstream file;
      fileMaster().openInputFile(filename, file);
      readFieldRGrid(file, field, unitCell);
      file.close();
   }

   template <int D>
   void FieldIo<D>::readFieldsRGridData(std::istream& in,
                                        DArray< RField<D> >& fields,
                                        int nMonomer)
   const
   {
      checkAllocateFields(fields, mesh().dimensions(), nMonomer);
      readRGridData(in, fields, mesh().dimensions(), nMonomer);
   }

   template <int D>
   void FieldIo<D>::writeFieldsRGrid(std::ostream &out,
                                     DArray<RField<D> > const & fields,
                                     UnitCell<D> const & unitCell,
                                     bool writeHeader,
                                     bool isSymmetric,
                                     bool writeMeshSize) const
   {
      int nMonomer = fields.capacity();
      UTIL_CHECK(nMonomer > 0);

      // Header
      if (writeHeader){
         writeFieldHeader(out, nMonomer, unitCell, isSymmetric);
      }
      if (writeMeshSize){
         IntVec<D> meshDimensions = fields[0].meshDimensions();
         writeMeshDimensions(out, meshDimensions);
      }

      // Data
      writeRGridData(out, fields, mesh().dimensions(), nMonomer);
   }

   template <int D>
   void FieldIo<D>::writeFieldsRGrid(std::string filename,
                                     DArray< RField<D> > const & fields,
                                     UnitCell<D> const & unitCell,
                                     bool isSymmetric) const
   {
      std::ofstream file;
      fileMaster().openOutputFile(filename, file);
      bool writeHeader = true;
      bool writeMeshSize = true;
      writeFieldsRGrid(file, fields, unitCell,
                       writeHeader, isSymmetric, writeMeshSize);
      file.close();
   }

   template <int D>
   void FieldIo<D>::writeFieldRGrid(std::ostream &out,
                                    RField<D> const & field,
                                    UnitCell<D> const & unitCell,
                                    bool writeHeader,
                                    bool isSymmetric) const
   {
      if (writeHeader) {
         writeFieldHeader(out, 1, unitCell, isSymmetric);
         writeMeshDimensions(out, mesh().dimensions());
      }
      writeRGridData(out, field, mesh().dimensions());
   }

   template <int D>
   void FieldIo<D>::writeFieldRGrid(std::string filename,
                                    RField<D> const & field,
                                    UnitCell<D> const & unitCell,
                                    bool isSymmetric) const
   {
      std::ofstream file;
      fileMaster().openOutputFile(filename, file);
      writeFieldRGrid(file, field, unitCell, isSymmetric);
      file.close();
   }

   // K-Grid Field Format

   template <int D>
   void FieldIo<D>::readFieldsKGrid(std::istream &in,
                                    DArray<RFieldDft<D> >& fields,
                                    UnitCell<D>& unitCell)
   const
   {
      int nMonomer;
      bool isSymmetric;
      FieldIo<D>::readFieldHeader(in, nMonomer, unitCell, isSymmetric);
      readMeshDimensions(in);
     
      checkAllocateFields(fields, mesh().dimensions(), nMonomer);
      readKGridData(in, fields, fields[0].dftDimensions(), nMonomer);
   }

   template <int D>
   void FieldIo<D>::readFieldsKGrid(std::string filename,
                                    DArray< RFieldDft<D> >& fields,
                                    UnitCell<D>& unitCell)
   const
   {
      std::ifstream file;
      fileMaster().openInputFile(filename, file);
      readFieldsKGrid(file, fields, unitCell);
      file.close();
   }

   template <int D>
   void FieldIo<D>::writeFieldsKGrid(std::ostream &out,
                                     DArray<RFieldDft<D> > const & fields,
                                     UnitCell<D> const & unitCell,
                                     bool isSymmetric) const
   {
      // Inspect fields array
      int nMonomer = fields.capacity();
      UTIL_CHECK(nMonomer > 0);
      for (int i = 0; i < nMonomer; ++i) {
         UTIL_CHECK(fields[i].meshDimensions() == mesh().dimensions());
      }
      IntVec<D> dftDimensions = fields[0].dftDimensions();

      // Write header
      writeFieldHeader(out, nMonomer, unitCell, isSymmetric);
      writeMeshDimensions(out, mesh().dimensions());

      writeKGridData(out, fields, dftDimensions, nMonomer);
   }

   template <int D>
   void FieldIo<D>::writeFieldsKGrid(std::string filename,
                                    DArray< RFieldDft<D> > const & fields,
                                    UnitCell<D> const & unitCell,
                                    bool isSymmetric) const
   {
      std::ofstream file;
      fileMaster().openOutputFile(filename, file);
      writeFieldsKGrid(file, fields, unitCell, isSymmetric);
      file.close();
   }

   // File IO Utilities

   /*
   * Read common part of field header and extract
   * the number of monomers (number of fields) in the file.
   */
   template <int D>
   void FieldIo<D>::readFieldHeader(std::istream& in,
                                    int& nMonomer,
                                    UnitCell<D>& unitCell,
                                    bool & isSymmetric) const
   {
      // Preconditions
      UTIL_CHECK(latticePtr_);
      if (unitCell.lattice() == UnitCell<D>::Null) {
         UTIL_CHECK(unitCell.nParameter() == 0);
      } else {
         UTIL_CHECK(unitCell.nParameter() > 0);
         UTIL_CHECK(unitCell.lattice() == lattice());
      }

      // Read field header to set unitCell, groupNameIn, nMonomer
      int ver1, ver2;
      std::string groupNameIn;

      Pscf::Prdc::readFieldHeader(in, ver1, ver2, unitCell,
                                  groupNameIn, nMonomer);
      // Note: Function definition in prdc/crystal/fieldHeader.tpp

      // Checks of data from header
      UTIL_CHECK(ver1 == 1);
      //UTIL_CHECK(ver2 == 0);
      UTIL_CHECK(unitCell.isInitialized());
      UTIL_CHECK(unitCell.lattice() != UnitCell<D>::Null);
      UTIL_CHECK(unitCell.nParameter() > 0);

      // Validate or initialize lattice type
      if (lattice() != unitCell.lattice()) {
         Log::file() << std::endl
               << "Error - Mismatched lattice types "
               << "in function FieldIo<D>::readFieldHeader:\n"
               << "  FieldIo::lattice  :" << lattice() << "\n"
               << "  Unit cell lattice :" << unitCell.lattice()
               << "\n";
         UTIL_THROW("Mismatched lattice types");
      }

      // Check for presence of group name
      isSymmetric = false;
      if (groupNameIn != "") {
         isSymmetric = true;
      }

      // Process group and basis (if any)
      if (hasGroup()) {

         // Check consistency of groupName values
         if (isSymmetric) {
            UTIL_CHECK(groupNamePtr_);
            if (groupNameIn != groupName()) {
               Log::file() << std::endl
                  << "Error - Mismatched group names in "
                  << "function FieldIo<D>::readFieldHeader:\n"
                  << "  FieldIo::groupName :" << groupName() << "\n"
                  << "  Field file header  :" << groupNameIn << "\n";
               UTIL_THROW("Mismatched group names");
            }
         }

         // If there is a group but no basis, construct a basis
         UTIL_CHECK(basisPtr_);
         if (!basis().isInitialized()) {
            basisPtr_->makeBasis(mesh(), unitCell, group());
         }
         UTIL_CHECK(basis().isInitialized());

      }

   }

   template <int D>
   void FieldIo<D>::readMeshDimensions(std::istream& in) const
   {
      // Read and check input stream mesh dimensions
      std::string label;
      in >> label;
      if (label != "mesh" && label != "ngrid") {
         std::string msg =  "\n";
         msg += "Error reading field file:\n";
         msg += "Expected mesh or ngrid, but found [";
         msg += label;
         msg += "]";
         UTIL_THROW(msg.c_str());
      }
      IntVec<D> nGrid;
      in >> nGrid;
      if (nGrid != mesh().dimensions()) {
         Log::file()
             << "Inconsistent mesh dimensions:\n"
             << "mesh().dimensions() = " << mesh().dimensions() << "\n"
             << "nGrid (from file)   = " << mesh().dimensions() << "\n";
         UTIL_THROW("Inconsistent dimensions in Mesh and field header");
      }
   }

   template <int D>
   void FieldIo<D>::writeFieldHeader(std::ostream &out,
                                     int nMonomer,
                                     UnitCell<D> const & unitCell,
                                     bool isSymmetric) const
   {
      int v1 = 1;
      int v2 = 0;
      std::string gName = "";
      if (isSymmetric) {
         UTIL_CHECK(hasGroup());
         gName = groupName();
      }
      Pscf::Prdc::writeFieldHeader(out, v1, v2, unitCell,
                                   gName, nMonomer);
      // Note: This function is defined in prdc/crystal/fieldHeader.tpp
   }

   template <int D>
   void
   FieldIo<D>::writeMeshDimensions(std::ostream &out,
                                   IntVec<D> const& meshDimensions)
   const
   {
      out << "mesh " <<  std::endl
          << "           " << meshDimensions << std::endl;
   }

   // Field Format Conversion Functions

   template <int D>
   void FieldIo<D>::convertBasisToKGrid(DArray<double> const & in,
                                        RFieldDft<D>& out) const
   {  Prdc::convertBasisToKGrid(in, out, basis(), out.dftDimensions()); }

   template <int D>
   void
   FieldIo<D>::convertBasisToKGrid(DArray< DArray <double> > const & in,
                                   DArray< RFieldDft<D> >& out) const
   {
      UTIL_ASSERT(in.capacity() == out.capacity());
      int n = in.capacity();
      for (int i = 0; i < n; ++i) {
         convertBasisToKGrid(in[i], out[i]);
      }
   }

   template <int D>
   void FieldIo<D>::convertKGridToBasis(RFieldDft<D> const & in,
                                        DArray<double>& out,
                                        bool checkSymmetry,
                                        double epsilon) const
   {
      Prdc::convertKGridToBasis(in, out, basis(), in.dftDimensions(),
                                checkSymmetry, epsilon);
   }

   template <int D>
   void 
   FieldIo<D>::convertKGridToBasis(DArray< RFieldDft<D> > const & in,
                                   DArray< DArray <double> > & out,
                                   bool checkSymmetry,
                                   double epsilon) const
   {
      UTIL_ASSERT(in.capacity() == out.capacity());
      int n = in.capacity();

      bool symmetric(true);
      for (int i = 0; i < n; ++i) {
         if (checkSymmetry) {
            // Check if kgrid has symmetry
            bool tmp_sym = hasSymmetry(in[i], epsilon, true);
            if (!tmp_sym) symmetric = false;
         }
         convertKGridToBasis(in[i], out[i], false);
      }

      // Print warning if any input field is not symmetric
      if (!symmetric) {
         Log::file() << std::endl
            << "WARNING: non-negligible error in conversion to "
            << "symmetry-adapted basis format." << std::endl
            << "See error values printed above for each asymmetric field."
            << std::endl
            << "The field that is output by this operation will be "
            << "a symmetrized version of" << std::endl
            << "the input field." << std::endl << std::endl;
      }

   }

   template <int D>
   void
   FieldIo<D>::convertBasisToRGrid(DArray<double> const & in,
                                   RField<D>& out) const
   {
      checkAllocateField(workDft_, mesh().dimensions());
      convertBasisToKGrid(in, workDft_);
      fft().inverseTransformSafe(workDft_, out);
   }

   template <int D>
   void
   FieldIo<D>::convertBasisToRGrid(DArray< DArray <double> > const & in,
                                   DArray< RField<D> >& out) const
   {
      UTIL_ASSERT(in.capacity() == out.capacity());
      checkAllocateField(workDft_, mesh().dimensions());

      int n = in.capacity();
      for (int i = 0; i < n; ++i) {
         convertBasisToKGrid(in[i], workDft_);
         fft().inverseTransformSafe(workDft_, out[i]);
      }
   }

   template <int D>
   void
   FieldIo<D>::convertRGridToBasis(RField<D> const & in,
                                   DArray<double> & out,
                                   bool checkSymmetry,
                                   double epsilon) const
   {
      checkAllocateField(workDft_, mesh().dimensions());
      fft().forwardTransform(in, workDft_);
      convertKGridToBasis(workDft_, out, checkSymmetry, epsilon);
   }

   template <int D>
   void
   FieldIo<D>::convertRGridToBasis(DArray< RField<D> > const & in,
                                   DArray< DArray <double> > & out,
                                   bool checkSymmetry,
                                   double epsilon) const
   {
      UTIL_ASSERT(in.capacity() == out.capacity());
      checkAllocateField(workDft_, mesh().dimensions());

      int n = in.capacity();

      bool symmetric(true);
      for (int i = 0; i < n; ++i) {
         fft().forwardTransform(in[i], workDft_);
         if (checkSymmetry) {
            // Check if kgrid has symmetry
            bool tmp_sym = hasSymmetry(workDft_, epsilon, true);
            if (!tmp_sym) symmetric = false;
         }
         convertKGridToBasis(workDft_, out[i], false);
      }

      // Print warning if any input fields is asymmetric
      if (!symmetric) {
         Log::file() << std::endl
             << "WARNING: non-negligible error in conversion to "
             << "symmetry-adapted basis format." << std::endl
             << "   See error values printed above for each "
             << "asymmetric field." << std::endl
             << "   The field that is output by the above operation "
             << "will be a" << std::endl
             << "   symmetrized version of the input field."
             << std::endl << std::endl;
      }
   }

   /*
   * Apply inverse FFT to an array of k-grid fields.
   */
   template <int D>
   void
   FieldIo<D>::convertKGridToRGrid(DArray< RFieldDft<D> > & in,
                                   DArray< RField<D> >& out) const
   {
      UTIL_ASSERT(in.capacity() == out.capacity());
      int n = in.capacity();
      for (int i = 0; i < n; ++i) {
         fft().inverseTransformSafe(in[i], out[i]);
      }
   }

   /*
   * Apply inverse FFT to a single k-grid field.
   */
   template <int D>
   void
   FieldIo<D>::convertKGridToRGrid(RFieldDft<D>& in, RField<D>& out) const
   {
      fft().inverseTransformSafe(in, out);
   }

   /*
   * Apply forward FFT to an array of r-grid fields.
   */
   template <int D>
   void
   FieldIo<D>::convertRGridToKGrid(DArray< RField<D> > const & in,
                                   DArray< RFieldDft<D> >& out) const
   {
      UTIL_ASSERT(in.capacity() == out.capacity());
      int n = in.capacity();
      for (int i = 0; i < n; ++i) {
         fft().forwardTransform(in[i], out[i]);
      }
   }

   /*
   * Apply forward FFT to a single r-grid field.
   */
   template <int D>
   void
   FieldIo<D>::convertRGridToKGrid(RField<D> const & in,
                                   RFieldDft<D>& out) const
   {
      fft().forwardTransform(in, out);
   }

   /*
   * Test if an RField<D> has declared space group symmetry.
   * Return true if symmetric, false otherwise. Print error values
   * if verbose == true and hasSymmetry == false.
   */
   template <int D>
   bool FieldIo<D>::hasSymmetry(RField<D> const & in, double epsilon,
                                bool verbose) const
   {
      checkAllocateField(workDft_, mesh().dimensions());
      fft().forwardTransform(in, workDft_);
      return hasSymmetry(workDft_, epsilon, verbose);
   }

   /*
   * Test if an RFieldDft has the declared space group symmetry.
   * Return true if symmetric, false otherwise. Print error values
   * if verbose == true and hasSymmetry == false.
   */
   template <int D>
   bool FieldIo<D>::hasSymmetry(RFieldDft<D> const & in, double epsilon,
                                bool verbose) const
   {
      return Prdc::hasSymmetry(in, basis(), in.dftDimensions(), 
                               epsilon, verbose);
   }

   // Grid Manipulation Utilities

   // Note: explicit instantiations of expandRGridDimension
   // are defined in FieldIo.cpp

   template <int D>
   void
   FieldIo<D>::expandRGridDimension(std::string filename,
                                    DArray< RField<D> > const & fields,
                                    UnitCell<D> const & unitCell,
                                    int d,
                                    DArray<int> newGridDimensions) const
   {
      std::ofstream file;
      fileMaster().openOutputFile(filename, file);
      expandRGridDimension(file, fields, unitCell, d, newGridDimensions);
      file.close();
   }

   template <int D>
   void FieldIo<D>::replicateUnitCell(std::ostream &out,
                                      DArray< RField<D> > const & fields,
                                      UnitCell<D> const & unitCell,
                                      IntVec<D> const & replicas) const

   {
      // Obtain number of monomer types
      int nMonomer = fields.capacity();
      UTIL_CHECK(nMonomer > 0);

      // Obtain initial dimensions of fields
      IntVec<D> meshDimensions = fields[0].meshDimensions();

      // Define dimension of replicated fields
      IntVec<D> replicateDimensions;
      for (int i = 0; i < D; ++i) {
         UTIL_CHECK(replicas[i] != 0);
         replicateDimensions[i] = replicas[i] * meshDimensions[i];
      }

      // Set up new UnitCell
      UnitCell<D> cell;
      FSArray<double, 6> parameters;
      int nParameter = unitCell.nParameter();
      for (int i = 0; i < nParameter; i++) {
         parameters[i]=  replicas[i]* unitCell.parameter(i);
      }
      cell.set(unitCell.lattice(), parameters);

      // Allocate outFields
      DArray<RField<D> > outFields;
      outFields.allocate(nMonomer);
      for (int i = 0; i < nMonomer; ++i) {
         outFields[i].allocate(replicateDimensions);
      }

      int n1 = 0;
      int n2 = 0;
      int n3 = 0;
      int ybeginPtr = 0;
      int zbeginPtr = 0;
      int rank = 0;

      if (D == 1) {
         parameters[0] = replicas[0] * unitCell.parameter(0);
         cell.set(unitCell.lattice(), parameters);
         for (int counter = 0; counter < replicas[0]; counter++) {
            for (int j = 0; j < meshDimensions[0]; j++) {
               for (int i = 0; i < nMonomer; i++) {
                  outFields[i][rank] = fields[i][j];
               }
               rank++;
            }
         }
      }

      if (D == 2) {
         for (int yCounter = 0; yCounter < replicas[0]; yCounter++) {
            n2 = 0;
            while (n1 < meshDimensions[0]) {
               ybeginPtr = n1 * meshDimensions[1];
               for (int xCounter = 0; xCounter < replicas[1]; xCounter++) {
                  n1 = 0;
                  int r = ybeginPtr;
                  while (n2 < meshDimensions[1]) {
                     for (int i = 0; i < nMonomer; i++) {
                        outFields[i][rank] = fields[i][r];
                     }
                     rank++;
                     n2++;
                     r++;
                  }
               }
               n1++;
            }
         }
      }

      if (D == 3) {
         for (int zCounter = 0; zCounter < replicas[0]; zCounter++) {
            n1 = 0;
            while (n1 < meshDimensions[0]) {
               zbeginPtr =  n1 * meshDimensions[2] * meshDimensions[1];
               for (int yCounter = 0; yCounter < replicas[1]; yCounter++) {
                  n2 = 0;
                  while (n2 < meshDimensions[1]) {
                     ybeginPtr = zbeginPtr + n2 * meshDimensions[2];
                     for (int xCounter = 0; xCounter < replicas[2]; xCounter++)
                     {
                        n3 = 0;
                        int r = ybeginPtr;
                        while (n3 < meshDimensions[2]) {
                           for (int i = 0; i < nMonomer; i++) {
                              outFields[i][rank] = fields[i][r];
                           }
                           rank++;
                           n3++;
                           r++;
                        }
                     }
                     n2++;
                  }
               }
               n1++;
            }
         }
      }

      // Write header
      int v1 = 1;
      int v2 = 0;
      std::string gName = "";
      Pscf::Prdc::writeFieldHeader(out, v1, v2, cell, gName, nMonomer);
      writeMeshDimensions(out, replicateDimensions);

      // Write field data
      writeRGridData(out, outFields, replicateDimensions, nMonomer);
   }

   template <int D>
   void FieldIo<D>::replicateUnitCell(std::string filename,
                                      DArray<RField<D> > const & fields,
                                      UnitCell<D> const & unitCell,
                                      IntVec<D> const & replicas) const
   {
      std::ofstream file;
      fileMaster().openOutputFile(filename, file);
      replicateUnitCell(file, fields, unitCell, replicas);
      file.close();
   }

} // namespace Rpc
} // namespace Pscf
#endif
