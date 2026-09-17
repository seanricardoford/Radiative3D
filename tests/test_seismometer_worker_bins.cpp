#include <cassert>
#include <iostream>

#include "dataout.hpp"
#include "ecs.hpp"
#include "phonons.hpp"
#include "simulation_report.hpp"

int main() {
  Seismometer::SetTimeBinParameters(1.0, 4.0);

  Real inner[RAY_NBT] = {0.5, 0.5};
  Real outer[RAY_NBT] = {2.0, 2.0};
  Seismometer seis(R3::XYZ(0.0, 0.0, 0.0),
                   R3::XYZ(1.0, 0.0, 0.0),
                   inner,
                   outer,
                   "ENZ");

  SimulationReportContext context(1, 4);
  Phonon phon(R3::XYZ(1.0, 0.0, 0.0),
              S2::ThetaPhi(Geometry::Pi * 0.5, 0.0),
              RAY_P);

  bool passthrough = seis.CatchPhonon(phon, context.BinsFor(0));

  assert(!passthrough);
  assert(context.BinsFor(0)[0].mCountByType[RAY_P] == 1);
  assert(context.BinsFor(0)[0].mEnergyByType[RAY_P] > 0.0);
  assert(context.BinsFor(0)[1].mCountByType[RAY_P] == 0);

  std::cout << "seismometer worker-bin tests passed\n";
  return 0;
}
