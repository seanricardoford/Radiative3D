#include <cassert>
#include <iostream>

#include "events.hpp"
#include "phonons.hpp"
#include "simulation_report.hpp"

int main() {
  typedef Phonon (ShearDislocation::*EventGenerator)(
      RandomEngine &, SimulationReportContext &);
  typedef void (Phonon::*Propagator)(RandomEngine &,
                                     SimulationReportContext &);

  EventGenerator generate = &ShearDislocation::GenerateEventPhonon;
  Propagator propagate = &Phonon::Propagate;

  assert(generate != 0);
  assert(propagate != 0);
  std::cout << "parallel context API tests passed\n";
  return 0;
}
