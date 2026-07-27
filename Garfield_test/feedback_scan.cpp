#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>

#include "Garfield/ComponentTcad2d.hh"
#include "Garfield/MediumSilicon.hh"
#include "Garfield/Sensor.hh"

#include "lgad_tools.hh"

using namespace Garfield;

namespace {

void PrintUsage(const char* exe) {
  std::cerr
      << "Usage:\n  " << exe
      << " FIELD.sta BIAS_V X_UM OUT_DIR [MODEL] [MIN_EVENTS]"
         " [MAX_EVENTS] [TARGET_REL_SEM] [FINE_STEP_NM] [CAP]\n\n"
      << "Example (deep-JTE 180 V):\n  " << exe
      << " /home/ahaines561/HEP/MAS/Silvaco_dat/lgad180V.sta"
         " 180 80 output_feedback_180V/ okuto 500 5000 0.05 50 5000\n";
}

}  // namespace

int main(int argc, char* argv[]) {
  if (argc < 5) {
    PrintUsage(argv[0]);
    return 2;
  }

  const std::string fieldFile = argv[1];
  const double biasV = std::atof(argv[2]);
  const double requestedXUm = std::atof(argv[3]);
  std::string outDir = argv[4];
  const std::string model = argc > 5 ? argv[5] : "okuto";
  const int minEvents = argc > 6 ? std::atoi(argv[6]) : 500;
  const int maxEvents = argc > 7 ? std::atoi(argv[7]) : 5000;
  const double targetRelSem = argc > 8 ? std::atof(argv[8]) : 0.05;
  const double fineStepNm = argc > 9 ? std::atof(argv[9]) : 50.;
  const std::size_t sizeCap = argc > 10
      ? static_cast<std::size_t>(std::strtoull(argv[10], nullptr, 10))
      : 5000;

  if (outDir.empty()) outDir = "./";
  if (outDir.back() != '/') outDir.push_back('/');
  std::filesystem::create_directories(outDir);

  std::cout << "field map = " << fieldFile << "\n"
            << "bias label = " << biasV << " V\n"
            << "requested x = " << requestedXUm << " um\n"
            << "model = " << model << "\n"
            << "output = " << outDir << "\n";

  ComponentTcad2d cmp;
  if (!cmp.InitialiseSilvaco(fieldFile)) return 1;

  MediumSilicon si;
  si.SetTemperature(293.15);
  if (!SetImpactIonisationModel(si, model)) return 2;
  DumpTownsendCoefficients(
      si, outDir + "townsend_garfield_" + model + "_293K.csv");

  // The uploaded structures use Silvaco region 1 as the silicon body,
  // corresponding to Garfield region index 0. Oxide and metal regions must
  // not be assigned a drift medium.
  if (cmp.GetNumberOfRegions() == 0) {
    std::cerr << "Field map contains no regions.\n";
    return 1;
  }
  cmp.SetMedium(0, &si);
  cmp.SetRangeZ(-5.e-4, 5.e-4);
  cmp.PrintRegions();

  double bx0 = 0., by0 = 0., bz0 = 0.;
  double bx1 = 0., by1 = 0., bz1 = 0.;
  if (!cmp.GetBoundingBox(bx0, by0, bz0, bx1, by1, bz1)) {
    std::cerr << "Could not get field-map bounding box.\n";
    return 1;
  }

  // Preserve the small mesh-boundary nudge used by the existing toolkit.
  double x0 = requestedXUm * 1.e-4 + 0.13e-4;
  if (x0 <= bx0 || x0 >= bx1) {
    std::cerr << "Requested x lies outside the field map: x range is ["
              << bx0 * 1.e4 << ", " << bx1 * 1.e4 << "] um.\n";
    return 2;
  }

  const std::string biasLabel = FormatBias(biasV);
  const auto profile = ScanFieldProfile(
      cmp, x0, by0, by1, 800,
      outDir + "profile_" + biasLabel + "V.csv", 100.);
  if (!profile.valid) return 1;

  std::cout << "actual probe x = " << x0 * 1.e4 << " um\n";

  Sensor sensor;
  sensor.AddComponent(&cmp);
  sensor.SetArea(bx0, profile.yTop + 0.02e-4, -5.e-4,
                 bx1, profile.yBot, 5.e-4);

  const double fineStepCm = fineStepNm * 1.e-7;
  const double bulkStepCm = 250. * 1.e-7;
  const double fineBandHalfWidthCm = 2.5e-4;
  const double yFineLo =
      std::max(profile.yGain - fineBandHalfWidthCm, profile.yTop);
  const double yFineHi =
      std::min(profile.yGain + fineBandHalfWidthCm, profile.yBot);

  // Inject above the gain peak, matching the existing ladder test. The
  // primary electron drifts through the multiplication region; the primary
  // hole moves away from it. This isolates secondary-hole feedback.
  const double yInj =
      std::min(profile.yGain + 5.e-4, 0.5 * (profile.yGain + profile.yBot));

  double ladderStepCm = fineStepCm;
  std::atomic<long long> calls{0}, fineCalls{0}, coarseCalls{0}, printed{0};
  AvalancheMC avalanche;
  avalanche.SetSensor(&sensor);
  ConfigureAvalanche(avalanche, ladderStepCm, bulkStepCm,
                     yFineLo, yFineHi, 6., sizeCap,
                     /*enableSignal=*/false,
                     /*enableMultithreading=*/false,
                     "feedback", calls, fineCalls, coarseCalls, printed);

  FeedbackScanConfig scan;
  scan.models = {model};
  scan.minEvents = minEvents;
  scan.maxEvents = maxEvents;
  scan.batchSize = std::max(50, std::min(250, minEvents));
  scan.targetRelativeSem = targetRelSem;
  scan.heavyTailThresholdF = 3.0;
  scan.runHoleSeed = true;
  scan.runPairSeed = true;

  RunModelComparison(avalanche, si, x0, yInj, ladderStepCm, sizeCap,
                     outDir, biasLabel, profile.eMax, profile.eMax, scan);

  std::cout << "[stepfn feedback] calls=" << calls.load()
            << " fine=" << fineCalls.load()
            << " coarse=" << coarseCalls.load() << "\n";
  return 0;
}
