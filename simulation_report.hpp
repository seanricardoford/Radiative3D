// simulation_report.hpp
//
// Worker-local simulation reporting state.  The propagation workers update
// these objects without sharing mutable counters or seismometer bins.  The
// owning DataReporter reduces them after all workers have joined.
//
#ifndef SIMULATION_REPORT_H_
#define SIMULATION_REPORT_H_

#include <cassert>
#include <cstddef>
#include <vector>

#include "raytype.hpp"
#include "typedefs.hpp"

struct SimulationReportBin {
  enum { NUM_ENERGY_AXES = 3 };

  Real mEnergyAxes[NUM_ENERGY_AXES];
  Real mEnergyByType[RAY_NUMBASICTYPES];
  unsigned mCountByType[RAY_NUMBASICTYPES];

  SimulationReportBin() {
    for (int axis = 0; axis < NUM_ENERGY_AXES; ++axis) {
      mEnergyAxes[axis] = 0;
    }
    for (int type = 0; type < RAY_NUMBASICTYPES; ++type) {
      mEnergyByType[type] = 0;
      mCountByType[type] = 0;
    }
  }

  void Add(const SimulationReportBin & other) {
    for (int axis = 0; axis < NUM_ENERGY_AXES; ++axis) {
      mEnergyAxes[axis] += other.mEnergyAxes[axis];
    }
    for (int type = 0; type < RAY_NUMBASICTYPES; ++type) {
      mEnergyByType[type] += other.mEnergyByType[type];
      mCountByType[type] += other.mCountByType[type];
    }
  }
};

class SimulationReportContext {
public:
  SimulationReportContext() :
    mSeismometerCount(0),
    mBinCount(0),
    mNumLost(0),
    mNumTimeout(0),
    mNumInvalid(0),
    mDiagInvalid(0) {}

  SimulationReportContext(std::size_t seismometer_count,
                          std::size_t bin_count) :
    mSeismometerCount(0),
    mBinCount(0),
    mNumLost(0),
    mNumTimeout(0),
    mNumInvalid(0),
    mDiagInvalid(0) {
    Reset(seismometer_count, bin_count);
  }

  void Reset(std::size_t seismometer_count, std::size_t bin_count) {
    mSeismometerCount = seismometer_count;
    mBinCount = bin_count;
    mBins.assign(mSeismometerCount * mBinCount, SimulationReportBin());
    mNumLost = 0;
    mNumTimeout = 0;
    mNumInvalid = 0;
    mDiagInvalid = 0;
  }

  SimulationReportBin * BinsFor(std::size_t seismometer_index) {
    assert(seismometer_index < mSeismometerCount);
    return &mBins[seismometer_index * mBinCount];
  }

  const SimulationReportBin * BinsFor(
      std::size_t seismometer_index) const {
    assert(seismometer_index < mSeismometerCount);
    return &mBins[seismometer_index * mBinCount];
  }

  void Merge(const SimulationReportContext & other) {
    assert(mSeismometerCount == other.mSeismometerCount);
    assert(mBinCount == other.mBinCount);
    assert(mBins.size() == other.mBins.size());
    for (std::size_t i = 0; i < mBins.size(); ++i) {
      mBins[i].Add(other.mBins[i]);
    }
    mNumLost += other.mNumLost;
    mNumTimeout += other.mNumTimeout;
    mNumInvalid += other.mNumInvalid;
    mDiagInvalid |= other.mDiagInvalid;
  }

  void RecordLost() { ++mNumLost; }
  void RecordTimeout() { ++mNumTimeout; }

  void RecordInvalid(unsigned reason) {
    assert(reason < sizeof(mDiagInvalid) * 8);
    ++mNumInvalid;
    mDiagInvalid |= (1u << reason);
  }

  std::size_t SeismometerCount() const { return mSeismometerCount; }
  std::size_t BinCount() const { return mBinCount; }
  unsigned long NumLost() const { return mNumLost; }
  unsigned long NumTimeout() const { return mNumTimeout; }
  unsigned long NumInvalid() const { return mNumInvalid; }
  unsigned InvalidDiagnostics() const { return mDiagInvalid; }

private:
  std::size_t mSeismometerCount;
  std::size_t mBinCount;
  std::vector<SimulationReportBin> mBins;
  unsigned long mNumLost;
  unsigned long mNumTimeout;
  unsigned long mNumInvalid;
  unsigned mDiagInvalid;
};

#endif  // SIMULATION_REPORT_H_
