#include <cassert>
#include <cstdint>
#include <iostream>

#include "probability.hpp"
#include "model.hpp"

int main() {
  const ModelParams defaults;
  assert(defaults.WorkerCount >= 1);
  assert(defaults.RandomSeed != 0);
  assert(defaults.ScatteringHorizontalLength == 0.0);
  assert(defaults.ScatteringVerticalLength == 0.0);

  const std::uint64_t seed = 0x123456789abcdef0ULL;

  RandomEngine first(seed);
  RandomEngine second(seed);
  for (int i = 0; i < 32; ++i) {
    assert(first.NextUInt64() == second.NextUInt64());
  }

  RandomEngine stream_a(RandomEngine::SeedForStream(seed, 7));
  RandomEngine stream_b(RandomEngine::SeedForStream(seed, 7));
  for (int i = 0; i < 32; ++i) {
    assert(stream_a.Uniform01() == stream_b.Uniform01());
  }

  RandomEngine stream_c(RandomEngine::SeedForStream(seed, 8));
  assert(stream_a.NextUInt64() != stream_c.NextUInt64());

  ProbDist distribution(3);
  distribution.SetRelativeProb(0, 1.0);
  distribution.SetRelativeProb(1, 2.0);
  distribution.SetRelativeProb(2, 3.0);
  RandomEngine chooser_a(seed);
  RandomEngine chooser_b(seed);
  for (int i = 0; i < 64; ++i) {
    assert(distribution.GetRandomIndex(chooser_a)
           == distribution.GetRandomIndex(chooser_b));
  }

  std::cout << "parallel feature tests passed\n";
  return 0;
}
