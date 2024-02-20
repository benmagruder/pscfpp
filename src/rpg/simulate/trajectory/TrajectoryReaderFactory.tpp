#ifndef RPG_TRAJECTORY_READER_TPP
#define RPG_TRAJECTORY_READER_TPP

#include "TrajectoryReaderFactory.h"

// Subclasses of ConfigIo
#include "FieldConfigReader.h"

namespace Pscf {
namespace Rpg {

   using namespace Util;

   /*
   * Constructor
   */
   template <int D>
   TrajectoryReaderFactory<D>::TrajectoryReaderFactory(System<D>& system)
    : sysPtr_(&system)
   {}

   /*
   * Return a pointer to a instance of Trajectory subclass className.
   */
   template <int D>
   TrajectoryReader<D>* 
   TrajectoryReaderFactory<D>::factory(const std::string &className) const
   {
      TrajectoryReader<D> *ptr = 0;

      // Try subfactories first
      ptr = trySubfactories(className);
      if (ptr) return ptr;

      if (className == "FieldConfigReader") {
        ptr = new FieldConfigReader<D>(*sysPtr_);
      } 
      return ptr;
   }

}
}
#endif
