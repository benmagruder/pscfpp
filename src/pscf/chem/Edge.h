#ifndef PSCF_EDGE_H
#define PSCF_EDGE_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "PolymerType.h"
#include <util/containers/Pair.h>

#include <iostream>

namespace Pscf
{

   using namespace Util;

   /**
   * Descriptor for a block within an acyclic block polymer.
   *
   * An Edge has:
   *
   *    - a monomer type id
   *    - a length
   *    - a block id (unique with the polymer)
   *    - ids of the two vertices at which it terminates
   *
   * An Edge is a descriptor for a block within a block polymer, but
   * does not provide functions or data structure needed to solve the
   * modified diffusion equation (MDE).  The Edge class is used as a 
   * base class for subclasses that are block MDE solvers as well as 
   * descriptors. Each implementation level sub-namespace of Pscf (i.e., 
   * R1d, Rpc and Rpg) contains a class named Block that is both a 
   * descriptor and an MDE solver for the associated block.  Each such 
   * Block class is a subclass of an instance of a class template
   *  Pscf::BlockTmpl, which is itself a subclass of the Edge class.
   *
   * Block objects associated with a polymer are normally stored in 
   * an array that is a private member of the Pscf::PolymerTmpl class 
   * template. The block id for each edge be set to the element index
   * of the Block within that array.  
   *
   * \ref user_param_block_sec "Parameter File Format"
   * \ingroup Pscf_Chem_Module
   */
   class Edge
   {
   public:

      /**
      * Constructor.
      */
      Edge();

      /**
      * Destructor.
      */
      virtual ~Edge();

      /**
      * Serialize to/from archive.
      *
      * \param ar  input or output Archive
      * \param versionId  archive format version index
      */
      template <class Archive>
      void serialize(Archive& ar, unsigned int versionId);

      /// \name Setters
      //@{

      /**
      * Set the id for this block.
      *
      * \param id  integer index for this block
      */
      void setId(int id);

      /**
      * Set indices of associated vertices.
      *
      * \param vertexAId  integer id of vertex A
      * \param vertexBId  integer id of vertex B
      */
      void setVertexIds(int vertexAId, int vertexBId);

      /**
      * Set the monomer type id.
      *
      * \param monomerId  integer id of monomer type
      */
      void setMonomerId(int monomerId);

      /**
      * Set the length of this block.
      *
      * The ``length" is steric volume / reference volume.
      *
      * \param length  block length (number of monomers).
      */
      virtual void setLength(double length);

      /**
      * Set the polymer type (branched or linear).
      *
      * By convention, if the polymer type of a block with block index id is
      * PolymerType::Linear, then vertexId(0) = id and vertexId(1) = id + 1.
      * The stream insertion and extraction operators for a Edge
      * can thus use a shorter string representation for linear polymers in
      * in which vertex ids are omitted.
      *
      * \param type  type of polymer (branched or linear)
      */
      void setPolymerType(PolymerType::Enum type);

      //@}
      /// \name Accessors (getters)
      //@{

      /**
      * Get the id of this block.
      */
      int id() const;

      /**
      * Get the monomer type id.
      */
      int monomerId() const;

      /**
      * Get the pair of associated vertex ids.
      */
      const Pair<int>& vertexIds() const;

      /**
      * Get id of an associated vertex.
      *
      * \param i index of vertex (0 or 1)
      */
      int vertexId(int i) const;

      /**
      * Get the length (number of monomers) in this block.
      */
      double length() const;

      /**
      * Get the type of the parent polymer (branched or linear).
      */
      PolymerType::Enum polymerType() const;

      //@}

   private:

      /// Identifier for this block, unique within the polymer.
      int id_;

      /// Identifier for the associated monomer type.
      int monomerId_;

      /// Length of this block = volume / monomer reference volume.
      double length_;

      /// Indexes of associated vertices
      Pair<int> vertexIds_;

      /// Type of polymer (branched or linear)
      PolymerType::Enum polymerType_;

      friend
      std::istream& operator >> (std::istream& in, Edge &block);

      friend
      std::ostream& operator << (std::ostream& out,
                                 const Edge &block);

   };

   /**
   * Input stream extractor (>>) for a Edge.
   *
   * Different text representations are used for linear and branched
   * polymers, based on the value of the polymerType enumeratin value.
   * Text representation for a branched polymer is:
   * \code
   *    monomerId length vertexId(0) vertexid(1)
   * \endcode
   * The vertex ids are omitted from the text representation for a linear
   * polymer, because the vertex ids for a linear polymer must be equal
   * to id and id + 1, where id denotes the block id.
   *
   * The polymerType must be set before a Edge can be read
   * from a stream. The block id must be set explicitly by calling
   * setId, rather than read from an istream. Vertex id values for blocks 
   * in a linear polymer must be set explicitly by calling setVertexIds
   * with consecutive values, as done in the function
   * Pscf::PolymerTmpl::readParameters.
   *
   * \param in  input stream
   * \param block  Edge to be read from stream
   * \return  modified input stream
   */
   std::istream& operator >> (std::istream& in, Edge &block);

   /**
   * Output stream inserter (<<) for a Edge.
   *
   * Different text representations are used for linear and branched
   * polymers, as discussed in documentation for the stream extractor
   * (>>) operator. Vertex ids are output only for blocks of branched
   * polymers.
   *
   * \param out  output stream
   * \param block  Edge to be written to stream
   * \return modified output stream
   */
   std::ostream&
   operator << (std::ostream& out, const Edge &block);

   // Inline member functions

   /*
   * Get the id of this block.
   */
   inline int Edge::id() const
   {  return id_; }

   /*
   * Get the monomer type id.
   */
   inline int Edge::monomerId() const
   {  return monomerId_; }

   /*
   * Get the pair of associated vertex ids.
   */
   inline const Pair<int>& Edge::vertexIds() const
   {  return vertexIds_; }

   /*
   * Get id of an associated vertex.
   */
   inline int Edge::vertexId(int i) const
   {  return vertexIds_[i]; }

   /*
   * Get the length (number of monomers) in this block.
   */
   inline double Edge::length() const
   {  return length_; }

   /*
   * Get the polymer type (branched or linear).
   */
   inline PolymerType::Enum Edge::polymerType() const
   {  return polymerType_; }

   /*
   * Serialize to/from an archive.
   */
   template <class Archive>
   void Edge::serialize(Archive& ar, unsigned int)
   {
      ar & id_;
      ar & monomerId_;
      ar & vertexIds_;
      ar & length_;
   }

}
#endif
