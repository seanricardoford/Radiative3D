#include <cassert>
#include <iostream>

#include "simulation_report.hpp"

int main() {
  SimulationReportContext first(2, 3);
  SimulationReportContext second(2, 3);

  SimulationReportBin * first_seis0 = first.BinsFor(0);
  first_seis0[1].mEnergyAxes[0] = 1.5;
  first_seis0[1].mEnergyByType[RAY_P] = 2.5;
  first_seis0[1].mCountByType[RAY_S] = 3;

  SimulationReportBin * second_seis0 = second.BinsFor(0);
  second_seis0[1].mEnergyAxes[0] = 4.5;
  second_seis0[1].mEnergyByType[RAY_P] = 5.5;
  second_seis0[1].mCountByType[RAY_S] = 7;
  second.BinsFor(1)[2].mEnergyAxes[2] = 8.5;

  first.RecordLost();
  first.RecordInvalid(2);
  second.RecordLost();
  second.RecordTimeout();
  second.RecordInvalid(5);

  first.Merge(second);

  assert(first.SeismometerCount() == 2);
  assert(first.BinCount() == 3);
  assert(first.NumLost() == 2);
  assert(first.NumTimeout() == 1);
  assert(first.NumInvalid() == 2);
  assert(first.InvalidDiagnostics() == ((1u << 2) | (1u << 5)));
  assert(first.BinsFor(0)[1].mEnergyAxes[0] == 6.0);
  assert(first.BinsFor(0)[1].mEnergyByType[RAY_P] == 8.0);
  assert(first.BinsFor(0)[1].mCountByType[RAY_S] == 10);
  assert(first.BinsFor(1)[2].mEnergyAxes[2] == 8.5);

  std::cout << "report reduction tests passed\n";
  return 0;
}
