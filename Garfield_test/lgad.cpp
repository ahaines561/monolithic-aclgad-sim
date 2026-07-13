#include <cmath>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <array>
#include <iterator>
#include <chrono>
#include <atomic>
#include <algorithm>
#include <iostream>
#include <TApplication.h>
#include <TCanvas.h>
#include <TSystem.h>

#include "Garfield/ComponentTcad2d.hh"
#include "Garfield/ComponentConstant.hh"
#include "Garfield/MediumSilicon.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/GeometrySimple.hh"
#include "Garfield/SolidBox.hh"
#include "Garfield/TrackHeed.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/ViewSignal.hh"
#include "Garfield/ViewDrift.hh"
#include "Garfield/ViewField.hh"

using namespace Garfield;
using Clock = std::chrono::steady_clock;

// diagnostic: seconds since t0, for per-stage wall-clock timing
double ElapsedS(const Clock::time_point& t0) {
  return std::chrono::duration<double>(Clock::now() - t0).count();
}

int main(int argc, char* argv[]) {
  const auto tRunStart = Clock::now();
  const bool doMIP = argc > 2 ? std::atoi(argv[2]) != 0 : true;
  // argv[3]: impact-ionisation model: vodm (default) | massey | grant | okuto
  const std::string model = argc > 3 ? argv[3] : "okuto";
  const double xTrackUm = argc > 4 ? std::atof(argv[4]) : 20.;
  // argv[5]: bias label for batch scans (e.g. "190"), written into
  // results.csv and per-model size-file names. Purely a label -- the
  // actual bias is whatever the .sta file was solved at.
  const std::string biasLabel = argc > 5 ? argv[5] : "NA";
  const std::string file = argc > 1 ? argv[1]
      : "/home/ahaines561/HEP/MAS/Silvaco_dat"
      // "/lgad.sta";
      // "/lgad180V.sta";
      "/lgad150V.sta";
  int rootArgc = 1;
  char* rootArgv[] = {argv[0], nullptr};
  TApplication app("app", &rootArgc, rootArgv);


  ComponentTcad2d cmp;
  if (!cmp.InitialiseSilvaco(file)) return 1;
  MediumSilicon si;
  // Canali model should be default, but set it explicitly
  si.SetSaturationVelocityModelCanali();
  si.SetDopingMobilityModelMasetti();
  si.SetTemperature(293.15);
  // gain is very sensitive to the ionisation model at these fields
  if (model == "massey") {
    si.SetImpactIonisationModelMassey();
  } else if (model == "grant") {
    si.SetImpactIonisationModelGrant();
  } else if (model == "okuto") {
    si.SetImpactIonisationModelOkutoCrowell();
  } else {
    si.SetImpactIonisationModelVanOverstraetenDeMan();
  }
  std::cout << "impact-ionisation model: "
            << (model == "massey" || model == "grant" || model == "okuto"
                    ? model : "vodm (van Overstraeten-de Man)")
            << std::endl;
  cmp.SetMedium("3", &si);
  cmp.SetRangeZ(-5.e-4, 5.e-4);
  cmp.PrintRegions();

  double bx0 = 0., by0 = 0., bz0 = 0., bx1 = 0., by1 = 0., bz1 = 0.;
  if (!cmp.GetBoundingBox(bx0, by0, bz0, bx1, by1, bz1)) {
    std::cerr << "Could not get the bounding box from the field map;\n"
              << "set the device extent manually in the macro.\n";
    return 1;
  }
  std::cout << "map extent: x = [" << bx0 * 1.e4 << ", " << bx1 * 1.e4
            << "] um, y = [" << by0 * 1.e4 << ", " << by1 * 1.e4
            << "] um" << std::endl;

  // dump dead points in the map (mid-silicon gap); compare against regions
  const auto tValidStart = Clock::now();
  {
    const int nx = 400, ny = 400;
    std::ofstream fvalid("validity_map.csv");
    fvalid << "x_um,y_um,status\n";
    std::size_t nInvalid = 0, nScanned = 0;
    for (int ix = 0; ix <= nx; ++ix) {
      const double x = bx0 + (bx1 - bx0) * ix / nx;
      for (int iy = 0; iy <= ny; ++iy) {
        const double y = by0 + (by1 - by0) * iy / ny;
        double ex, ey, ez, v; int st; Medium* m = nullptr;
        cmp.ElectricField(x, y, 0., ex, ey, ez, v, m, st);
        ++nScanned;
        if (st != 0) {
          ++nInvalid;
          fvalid << x * 1.e4 << "," << y * 1.e4 << "," << st << "\n";
        }
      }
    }
    std::cout << "validity scan: " << nInvalid << "/" << nScanned
              << " points invalid; written to validity_map.csv"
              << std::endl;
    std::cout << "  -> plot validity_map.csv (e.g. scatter x vs y) and "
              << "compare against the region geometry before deciding "
              << "whether a SetMedium() call is missing." << std::endl;
  }
  std::cout << "[timer] validity scan: " << ElapsedS(tValidStart)
            << " s" << std::endl;

  double x0 = xTrackUm * 1.e-4 + 0.13e-4;
  if (x0 <= bx0 || x0 >= bx1) {
    std::cout << "requested x = " << xTrackUm
              << " um is outside the map; using the device centre."
              << std::endl;
    x0 = 0.5 * (bx0 + bx1) + 0.13e-4;
  }
  std::cout << "probe/track line at x = " << x0 * 1.e4 << " um"
            << std::endl;

  std::cout << "# y[um]   V[V]   Ey[V/cm]\n";
  double eMax = 0., yGain = 0.;
  double yTop = 1., yBot = -1.;
  const int nScan = 400;
  for (int i = 0; i <= nScan; ++i) {
    const double y = by0 + ((by1 - by0) * i) / nScan;
    double ex, ey, ez, v; int st; Medium* m = nullptr;
    cmp.ElectricField(x0, y, 0., ex, ey, ez, v, m, st);
    if (st != 0) continue; 
    if (yTop > yBot) yTop = y;
    yBot = y;
    const double e = std::sqrt(ex * ex + ey * ey);
    if (e > eMax) { eMax = e; yGain = y; }
    std::cout << y * 1.e4 << "  " << v << "  " << ey << "\n";
  }
  if (yTop > yBot) {
    std::cerr << "No valid drift medium found along the scan line --\n"
              << "check the region/material assignment above.\n";
    return 1;
  }
  const double d = yBot - yTop;
  std::cout << "active silicon: y = [" << yTop * 1.e4 << ", "
            << yBot * 1.e4 << "] um (d = " << d * 1.e4 << " um)"
            << std::endl;
  std::cout << "Peak field " << eMax << " V/cm at y = " << yGain * 1.e4
            << " um  (gain layer)" << std::endl;
  if (eMax < 2.5e5) {
    std::cout << "NOTE: peak field < 250 kV/cm -- expect gain near 1 "
              << "at this bias." << std::endl;
  }

  double gMax = eMax, gx = x0, gy = yGain;
  const auto tHotspotStart = Clock::now();
  for (int ix = 0; ix <= 60; ++ix) {
    const double x = bx0 + ((bx1 - bx0) * ix) / 60. + 0.077e-4;
    for (int iy = 0; iy <= 150; ++iy) {
      const double y = by0 + ((by1 - by0) * iy) / 150. + 0.0113e-4;
      double ex, ey, ez, v; int st; Medium* m = nullptr;
      cmp.ElectricField(x, y, 0., ex, ey, ez, v, m, st);
      if (st != 0) continue;
      const double e = std::sqrt(ex * ex + ey * ey);
      if (e > gMax) { gMax = e; gx = x; gy = y; }
    }
    for (int iy = 0; iy <= 120; ++iy){
      const double y = yTop - 0.5e-4 + (3.5e-4 * iy) / 120.;
      double ex, ey, ez, v; int st; Medium* m = nullptr;
      cmp.ElectricField(x, y, 0., ex, ey, ez, v, m, st);
      if (st != 0) continue;
      const double e = std::sqrt(ex * ex + ey * ey);
      if (e > gMax) { gMax = e; gx = x; gy = y; }
    }
  }
  std::cout << "global max field " << gMax << " V/cm at (x, y) = ("
            << gx * 1.e4 << ", " << gy * 1.e4 << ") um" << std::endl;
  std::cout << "[timer] hotspot scan: " << ElapsedS(tHotspotStart)
            << " s" << std::endl;
  if (gMax > 1.2 * eMax) {
    std::cout << "NOTE: the global maximum exceeds the peak on the "
              << "configured line. To probe the hot spot instead, pass "
              << "its x on the command line, e.g.: ./lgad <file> <doMIP> "
              << "<model> " << gx * 1.e4 << std::endl;
  }

  // field dump over x = xTrack +- 5 um, full depth, for offline plotting
  {
    const auto tDumpStart = Clock::now();
    const double xLo = std::max(bx0, (xTrackUm - 5.) * 1.e-4);
    const double xHi = std::min(bx1, (xTrackUm + 5.) * 1.e-4);
    const double dxCm = 0.1e-4;
    const double dyCm = 0.05e-4; 
    const double eps = 0.013e-4;
    std::ofstream fdump("efield_window.csv");
    fdump << "x_um,y_um,Ex_Vcm,Ey_Vcm,Emag_Vcm,V\n";
    std::size_t nRows = 0, nBad = 0;
    double wMax = 0., wMaxX = 0., wMaxY = 0.;
    for (double x = xLo + eps; x <= xHi; x += dxCm) {
      for (double y = by0; y <= by1; y += dyCm) {
        double ex, ey, ez, v; int st; Medium* m = nullptr;
        cmp.ElectricField(x, y, 0., ex, ey, ez, v, m, st);
        if (st != 0) { ++nBad; continue; }
        const double e = std::sqrt(ex * ex + ey * ey);
        if (e > wMax) { wMax = e; wMaxX = x; wMaxY = y; }
        fdump << x * 1.e4 << "," << y * 1.e4 << "," << ex << "," << ey
              << "," << e << "," << v << "\n";
        ++nRows;
      }
    }
    std::cout << "field dump: x = [" << xLo * 1.e4 << ", " << xHi * 1.e4
              << "] um -> efield_window.csv (" << nRows << " points, "
              << nBad << " invalid skipped)" << std::endl;
    std::cout << "  window max |E| = " << wMax << " V/cm at (x, y) = ("
              << wMaxX * 1.e4 << ", " << wMaxY * 1.e4 << ") um"
              << std::endl;
    std::cout << "[timer] field dump: " << ElapsedS(tDumpStart) << " s"
              << std::endl;
  }

  {
    const auto tFullStart = Clock::now();
    const int nxFull = 250, nyFull = 250;
    const double eps = 0.013e-4;
    std::ofstream ffull("efield_full.csv");
    ffull << "x_um,y_um,Ex_Vcm,Ey_Vcm,Emag_Vcm,V\n";
    std::size_t nRows = 0, nBad = 0;
    for (int ix = 0; ix <= nxFull; ++ix) {
      const double x = bx0 + (bx1 - bx0) * ix / nxFull + eps;
      for (int iy = 0; iy <= nyFull; ++iy) {
        const double y = by0 + (by1 - by0) * iy / nyFull;
        double ex, ey, ez, v; int st; Medium* m = nullptr;
        cmp.ElectricField(x, y, 0., ex, ey, ez, v, m, st);
        if (st != 0) { ++nBad; continue; }
        const double e = std::sqrt(ex * ex + ey * ey);
        ffull << x * 1.e4 << "," << y * 1.e4 << "," << ex << "," << ey
              << "," << e << "," << v << "\n";
        ++nRows;
      }
    }
    std::cout << "full-device dump -> efield_full.csv (" << nRows
              << " points, " << nBad << " invalid skipped)" << std::endl;
    std::cout << "[timer] full-device dump: " << ElapsedS(tFullStart)
              << " s" << std::endl;
  }

  {
    const auto tHeedStart = Clock::now();
    SolidBox hbox(0.5 * (bx0 + bx1), 0.5 * (yTop + yBot), 0.,
                  0.5 * (bx1 - bx0), 0.5 * d, 5.e-4);
    GeometrySimple hgeo;
    hgeo.AddSolid(&hbox, &si);
    ComponentConstant hcmp;
    hcmp.SetGeometry(&hgeo);
    hcmp.SetElectricField(0., 100., 0.);
    Sensor hsensor;
    hsensor.AddComponent(&hcmp);
    TrackHeed htrack;
    htrack.SetSensor(&hsensor);
    htrack.SetParticle("pi");
    htrack.SetMomentum(180.e9);
    const int nHeedTracks = 300;
    std::ofstream fprim(biasLabel == "NA" ? "primary_pairs.txt"
        : ("primary_pairs_" + biasLabel + "V.txt"));
    for (int it = 0; it < nHeedTracks; ++it) {
      htrack.NewTrack(x0, yTop + 0.01e-4, 0., 0., 0., 1., 0.);
      double xc, yc, zc, tc, ec, extra;
      int nc = 0;
      unsigned long nPrimary = 0;
      while (htrack.GetCluster(xc, yc, zc, tc, nc, ec, extra)) {
        nPrimary += nc;
      }
      fprim << nPrimary << "\n";
    }
    std::cout << "wrote " << nHeedTracks << " fast Heed tracks (primary "
              << "pairs only, no avalanche) to primary_pairs"
              << (biasLabel == "NA" ? "" : "_" + biasLabel + "V")
              << ".txt" << std::endl;
    std::cout << "[timer] Heed primary-statistics: "
              << ElapsedS(tHeedStart) << " s" << std::endl;
  }

  // weighting field approx.
  ComponentConstant wcmp;
  wcmp.SetArea(bx0, yTop, -5.e-4, bx1, yBot, 5.e-4);
  wcmp.SetMedium(&si);
  wcmp.SetElectricField(0., 0., 0.);
  wcmp.SetWeightingField(0., 1. / d, 0., "pad");
  wcmp.SetWeightingPotential(0.5 * (bx0 + bx1), yTop, 0., 1.);
  {
    const double yMid = 0.5 * (yTop + yBot);
    double wx = 0., wy = 0., wz = 0.;
    wcmp.WeightingField(x0, yMid, 0., wx, wy, wz, "pad");
    std::cout << "Ew at centre = (" << wx << ", " << wy << ", " << wz
              << ")  [expect (0, " << 1. / d << ", 0)]" << std::endl;
    std::cout << "wpot: top = "
              << wcmp.WeightingPotential(x0, yTop, 0., "pad")
              << ", mid = "
              << wcmp.WeightingPotential(x0, yMid, 0., "pad")
              << ", back = "
              << wcmp.WeightingPotential(x0, yBot, 0., "pad")
              << "  [expect 1, 0.5, 0]" << std::endl;
  }

  Sensor sensor;
  sensor.AddComponent(&cmp);
  sensor.AddElectrode(&wcmp, "pad");
  sensor.SetTimeWindow(0., 0.005, 800);
  // extend to yBot -- was cutting the deepest (largest) clusters
  sensor.SetArea(bx0, yTop + 0.02e-4, -5.e-4, bx1, yBot, 5.e-4);

  // gain layer is ~0.5 um wide -- time steps under-resolved it, use distance steps
  // fine steps only pay off in the gain region; coarse in bulk (~50x faster)
  double kDistanceStepCm = 2.e-6;
  const double bulkStepCm = 2.5e-5;
  const double yFine = yGain + 2.5e-4;
  //verify the step-function args really are (x, y, z) in cm,
  std::atomic<long long> avalFnCalls{0}, avalFnFine{0}, avalFnCoarse{0},
      avalFnPrinted{0};
  AvalancheMC aval;
  aval.EnableMultithreading(true);
  aval.SetSensor(&sensor);
  aval.EnableSignalCalculation();
  aval.SetStepDistanceFunction([=, &avalFnCalls, &avalFnFine, &avalFnCoarse,
                                &avalFnPrinted](double x, double y,
                                                double z) {
    ++avalFnCalls;
    if (avalFnPrinted.fetch_add(1) < 5) {
      std::cout << "[stepfn aval] raw args: (" << x * 1.e4 << ", "
                << y * 1.e4 << ", " << z * 1.e4 << ") um  (interpreted as "
                << "x, y, z)" << std::endl;
    }
    if (y < yFine) { ++avalFnFine; return kDistanceStepCm; }
    ++avalFnCoarse;
    return bulkStepCm;
  });
  const std::size_t sizeCap = 5000;

  ViewDrift driftView;
  driftView.SetArea(bx0, yTop, -5.e-4, bx1, yBot, 5.e-4);
  aval.EnablePlotting(&driftView);

  // ViewField wfieldView;
  // wfieldView.SetComponent(&wcmp);
  // wfieldView.SetArea(0., -1.15308e-4, 0., 100.e-4, 50.e-4, 0.);

  const double yInj = std::min(yGain + 5.e-4, 0.5 * (yGain + yBot));
  std::atomic<long long> ladderFnCalls{0}, ladderFnFine{0}, ladderFnCoarse{0},
      ladderFnPrinted{0};
  AvalancheMC avalLadder;
  avalLadder.EnableMultithreading(true);
  avalLadder.SetSensor(&sensor);
  avalLadder.EnableSignalCalculation(false);  // ladder only needs sizes
  // divergent (f>=1) events dominate runtime -- cap the ladder sooner
  const std::size_t ladderCap = 5000;
  avalLadder.EnableAvalancheSizeLimit(ladderCap);
  aval.EnableAvalancheSizeLimit(sizeCap);
  // convergence scan of the gain-region step; bulk step stays at bulkStepCm
  double ladderStepCm = 2.e-6;  // 20 nm production step (from convergence scan)
  avalLadder.SetStepDistanceFunction(
      [&ladderStepCm, yFine, bulkStepCm, &ladderFnCalls, &ladderFnFine,
       &ladderFnCoarse, &ladderFnPrinted](double x, double y, double z) {
        ++ladderFnCalls;
        if (ladderFnPrinted.fetch_add(1) < 5) {
          std::cout << "[stepfn ladder] raw args: (" << x * 1.e4 << ", "
                    << y * 1.e4 << ", " << z * 1.e4 << ") um" << std::endl;
        }
        if (y < yFine) { ++ladderFnFine; return ladderStepCm; }
        ++ladderFnCoarse;
        return bulkStepCm;
      });
#if 0  // convergence scan (retired): G_eh plateau from 50 nm down; 20 nm adopted
  // convergence scan: pick the coarsest step where G stops changing
  const double distanceStepsCm[] = {1.e-5, 5.e-6, 2.e-6, 1.e-6, 5.e-7};
  std::vector<std::size_t> sizes;      // pooled only from the finest step
  std::size_t nCapped = 0;
  bool ehDivergent = false;
  std::cout << "# step[nm]   G_e (mean+-sem, N)      G_eh (mean+-sem, N)\n";
  const auto tLadderStart = Clock::now();
  for (std::size_t iStep = 0; iStep < std::size(distanceStepsCm); ++iStep) {
    const double stepCm = distanceStepsCm[iStep];
    const bool isFinest = (iStep + 1 == std::size(distanceStepsCm));
    ladderStepCm = stepCm;
    double res[2][2] = {{0., 0.}, {0., 0.}};   // [mode][sum, sum2]
    int nDone[2] = {0, 0};
    int nCapRung[2] = {0, 0};
    const int nWant[2] = {200, 300};           // high-stat pass
    for (int mode = 0; mode < 2; ++mode) {
      if (mode == 1 && ehDivergent) {
        std::cout << "step = " << stepCm * 1.e7
                  << " nm   G_eh skipped (divergent at coarser step)"
                  << std::endl;
        continue;
      }
      const auto tModeStart = Clock::now();
      for (int i = 0; i < nWant[mode]; ++i) {
        if (mode == 1 && i % 25 == 0) {
          std::cout << "  e+h injection " << i << std::endl;
        }
        const double xi = x0 + (i % 5 - 2) * 2.e-4;
        if (mode == 0) {
          avalLadder.AvalancheElectron(xi, yInj, 0., 0.);
        } else {
          avalLadder.AvalancheElectronHole(xi, yInj, 0., 0.);
        }
        std::size_t ne = 0, ni = 0;
        avalLadder.GetAvalancheSize(ne, ni);
        res[mode][0] += ne;
        res[mode][1] += double(ne) * double(ne);
        ++nDone[mode];
        if (ne >= ladderCap) ++nCapRung[mode];
        // only pool the finest step -- coarser ones aren't converged
        if (mode == 1 && isFinest) {
          sizes.push_back(ne);
          if (ne >= ladderCap) ++nCapped;
        }
      }
      // Print each mode's result as soon as it completes, so a
      // wedged avalanche in the other mode cannot hide it.
      const int n = nDone[mode];
      const double mean = res[mode][0] / n;
      const double var = res[mode][1] / n - mean * mean;
      const double sem = std::sqrt(var > 0. ? var / n : 0.);
      std::cout << "step = " << stepCm * 1.e7 << " nm   "
                << (mode == 0 ? "G_e  = " : "G_eh = ")
                << mean << " +- " << sem << " (N=" << n << ", capped="
                << nCapRung[mode] << ")"
                << "  [timer] " << ElapsedS(tModeStart) << " s"
                << std::endl;
      if (mode == 1 && nCapRung[1] * 10 > nWant[1]) {
        ehDivergent = true;
        std::cout << "NOTE: >10% of e+h avalanches hit the size cap ("
                  << ladderCap << ") -- hole-feedback divergence (f >= 1) "
                  << "at this bias for this ionisation model. G_eh is not "
                  << "a defined quantity here; skipping remaining e+h "
                  << "rungs. Compare models via G_e at this bias, or "
                  << "re-solve the Silvaco deck at lower bias."
                  << std::endl;
      }
    }
  }
  std::cout << "[stepfn ladder] calls=" << ladderFnCalls.load()
            << " fine=" << ladderFnFine.load() << " coarse="
            << ladderFnCoarse.load() << std::endl;
  std::cout << "[timer] full ladder: " << ElapsedS(tLadderStart)
            << " s" << std::endl;

  if (!sizes.empty()) {
    std::size_t mx = 0, nBig = 0;
    double mean = 0., mean2 = 0.;
    for (const auto s : sizes) { mean += s; mean2 += double(s) * double(s); }
    mean /= sizes.size();
    mean2 /= sizes.size();
    for (const auto s : sizes) {
      if (s > mx) mx = s;
      if (s > 5. * mean) ++nBig;
    }
    // excess noise factor F = <G^2>/<G>^2
    const double excessNoiseF = mean > 0. ? mean2 / (mean * mean) : 0.;
    std::cout << "e+h avalanche tail (finest step, " << sizes.size()
              << " events): max = " << mx << ", "
              << nBig << "/" << sizes.size() << " above 5x mean, "
              << nCapped << " capped at " << ladderCap
              << " (feedback-divergent)" << std::endl;
    std::cout << "excess noise factor F = <G^2>/<G>^2 = " << excessNoiseF
              << "  (McIntyre low-k limit -> F -> 1; large F/heavy tail "
              << "suggests operation near the breakdown knee)"
              << std::endl;
    std::ofstream fs("eh_sizes.txt");
    for (const auto s : sizes) fs << s << "\n";
    std::cout << "avalanche sizes (finest step only) written to "
              << "eh_sizes.txt" << std::endl;
  }
#endif

  // model comparison at the production step, same field map + injection points
  const char* cmpModels[] = {"vodm", "okuto", "massey", "grant"};
  const int nCmpModels = 4;
  struct ModelResult {
    double ge = 0., geSem = 0., geh = 0., gehSem = 0., F = 0., med = 0.;
    int nCapEh = 0, nEh = 0;
    bool divergent = false;
  } cmpRes[4];
  const auto tCmpStart = Clock::now();
  std::cout << "# model comparison: fine step = " << ladderStepCm * 1.e7
            << " nm, bulk step = " << bulkStepCm * 1.e7 << " nm, cap = "
            << ladderCap << std::endl;
  for (int im = 0; im < nCmpModels; ++im) {
    const std::string m = cmpModels[im];
    if (m == "massey") {
      si.SetImpactIonisationModelMassey();
    } else if (m == "okuto") {
      si.SetImpactIonisationModelOkutoCrowell();
    } else if (m == "grant") {
      si.SetImpactIonisationModelGrant();
    } else {
      si.SetImpactIonisationModelVanOverstraetenDeMan();
    }
    ModelResult& r = cmpRes[im];
    // G_e: single-pass electron gain (no hole transport, always finite)
    const auto tGeStart = Clock::now();
    double sum = 0., sum2 = 0.;
    const int nWantE = 200;
    for (int i = 0; i < nWantE; ++i) {
      const double xi = x0 + (i % 5 - 2) * 2.e-4;
      avalLadder.AvalancheElectron(xi, yInj, 0., 0.);
      std::size_t ne = 0, ni = 0;
      avalLadder.GetAvalancheSize(ne, ni);
      sum += ne;
      sum2 += double(ne) * double(ne);
    }
    r.ge = sum / nWantE;
    {
      const double var = sum2 / nWantE - r.ge * r.ge;
      r.geSem = std::sqrt(var > 0. ? var / nWantE : 0.);
    }
    std::cout << m << "   G_e  = " << r.ge << " +- " << r.geSem
              << " (N=" << nWantE << ")  [timer] "
              << ElapsedS(tGeStart) << " s" << std::endl;
    // G_eh: full gain with hole feedback; diverges if f >= 1
    const auto tEhStart = Clock::now();
    sum = 0.;
    sum2 = 0.;
    std::vector<std::size_t> mSizes;
    const int nWantEh = 300;
    for (int i = 0; i < nWantEh; ++i) {
      if (i % 50 == 0) std::cout << "  e+h injection " << i << std::endl;
      const double xi = x0 + (i % 5 - 2) * 2.e-4;
      avalLadder.AvalancheElectronHole(xi, yInj, 0., 0.);
      std::size_t ne = 0, ni = 0;
      avalLadder.GetAvalancheSize(ne, ni);
      sum += ne;
      sum2 += double(ne) * double(ne);
      mSizes.push_back(ne);
      if (ne >= ladderCap) ++r.nCapEh;
      // early exit: 50 events in, >10% capped -> feedback divergence
      if (i == 49 && r.nCapEh > 5) {
        r.divergent = true;
        std::cout << m << ": " << r.nCapEh << "/50 e+h avalanches hit the "
                  << "cap -- hole-feedback divergence (f >= 1) at this "
                  << "bias; G_eh undefined, aborting this model's e+h "
                  << "pass." << std::endl;
        break;
      }
    }
    r.nEh = int(mSizes.size());
    r.geh = sum / r.nEh;
    {
      const double var = sum2 / r.nEh - r.geh * r.geh;
      r.gehSem = std::sqrt(var > 0. ? var / r.nEh : 0.);
    }
    {
      double m1 = 0., m2 = 0.;
      for (const auto s : mSizes) { m1 += s; m2 += double(s) * double(s); }
      m1 /= r.nEh;
      m2 /= r.nEh;
      r.F = m1 > 0. ? m2 / (m1 * m1) : 0.;
      std::vector<std::size_t> tmp(mSizes);
      std::nth_element(tmp.begin(), tmp.begin() + tmp.size() / 2, tmp.end());
      r.med = double(tmp[tmp.size() / 2]);
    }
    if (!r.divergent && r.nCapEh * 10 > r.nEh) r.divergent = true;
    std::cout << m << "   G_eh = " << r.geh << " +- " << r.gehSem
              << " (N=" << r.nEh << ", median=" << r.med << ", capped="
              << r.nCapEh << (r.divergent ? ", DIVERGENT" : "")
              << ")   F = " << r.F << "  [timer] " << ElapsedS(tEhStart)
              << " s" << std::endl;
    std::ofstream fs(biasLabel == "NA" ? ("eh_sizes_" + m + ".txt")
                                        : ("eh_sizes_" + m + "_" + biasLabel
                                           + "V.txt"));
    for (const auto s : mSizes) fs << s << "\n";
  }
  std::cout << "[stepfn ladder] calls=" << ladderFnCalls.load()
            << " fine=" << ladderFnFine.load() << " coarse="
            << ladderFnCoarse.load() << std::endl;
  std::cout << "[timer] model comparison: " << ElapsedS(tCmpStart)
            << " s" << std::endl;
  std::cout << "\n# model   G_e            G_eh             median   F"
            << "      status" << std::endl;
  for (int im = 0; im < nCmpModels; ++im) {
    const ModelResult& r = cmpRes[im];
    std::cout << cmpModels[im] << "   " << r.ge << " +- " << r.geSem
              << "   " << r.geh << " +- " << r.gehSem << "   " << r.med
              << "   " << r.F << "   "
              << (r.divergent ? "DIVERGENT (G_eh unreliable)" : "ok")
              << std::endl;
  }

  // append one row per model to results.csv for the bias/model scan;
  // header written only if the file doesn't already exist
  {
    std::ifstream ftest("results.csv");
    const bool exists = ftest.good();
    ftest.close();
    std::ofstream fres("results.csv", std::ios::app);
    if (!exists) {
      fres << "bias,model,Epeak_line_Vcm,Epeak_global_Vcm,Ge,GeSem,Geh,"
           << "GehSem,median,F,nCapEh,N,divergent\n";
    }
    for (int im = 0; im < nCmpModels; ++im) {
      const ModelResult& r = cmpRes[im];
      fres << biasLabel << "," << cmpModels[im] << "," << eMax << ","
           << gMax << "," << r.ge << "," << r.geSem << "," << r.geh << ","
           << r.gehSem << "," << r.med << "," << r.F << "," << r.nCapEh
           << "," << r.nEh << "," << (r.divergent ? 1 : 0) << "\n";
    }
    std::cout << "appended " << nCmpModels << " rows to results.csv "
              << "(bias=" << biasLabel << ")" << std::endl;
  }

  // restore the CLI-selected model for the MIP section
  if (model == "massey") {
    si.SetImpactIonisationModelMassey();
  } else if (model == "grant") {
    si.SetImpactIonisationModelGrant();
  } else if (model == "okuto") {
    si.SetImpactIonisationModelOkutoCrowell();
  } else {
    si.SetImpactIonisationModelVanOverstraetenDeMan();
  }

  if (!doMIP) return 0;

  // MIP
  SolidBox box(0.5 * (bx0 + bx1), 0.5 * (yTop + yBot), 0.,
               0.5 * (bx1 - bx0), 0.5 * d, 5.e-4);
  GeometrySimple geo;
  geo.AddSolid(&box, &si);
  ComponentConstant cmpHeed;
  cmpHeed.SetGeometry(&geo);
  cmpHeed.SetElectricField(0., 100., 0.);
  Sensor heedSensor;
  heedSensor.AddComponent(&cmpHeed);

  TrackHeed track;
  track.SetSensor(&heedSensor);
  track.SetParticle("pi");
  track.SetMomentum(180.e9); 
  sensor.ClearSignal();

  track.NewTrack(x0, yTop + 0.01e-4, 0., 0., 0., 1., 0.); 
  double xc, yc, zc, tc, ec, extra;
  int nc = 0;
  unsigned long nPrimary = 0, nTotal = 0;
  unsigned long ncl = 0, nPairsDone = 0;
  double yLast = 0.;
  // save primaries to replay with multiplication off (no-gain reference)
  std::vector<std::array<double, 4>> primaries;  // {x, y, z, t}
  const auto tMipOnStart = Clock::now();
  while (track.GetCluster(xc, yc, zc, tc, nc, ec, extra)) {
    if (++ncl % 10 == 0) std::cout << "  cluster " << ncl << std::endl;
    yLast = yc;
    nPrimary += nc;
    for (int k = 0; k < nc; ++k) {
      primaries.push_back({xc, yc, zc, tc});
      aval.AvalancheElectronHole(xc, yc, zc, tc);
      std::size_t ne = 0, ni = 0;
      aval.GetAvalancheSize(ne, ni);
      nTotal += ne;
      // pair-level progress: clusters alone can hide long silent gaps,
      // since a single cluster's primary count varies widely (delta rays)
      if (++nPairsDone % 200 == 0) {
        std::cout << "  pair " << nPairsDone << "  ["
                  << ElapsedS(tMipOnStart) << " s elapsed]" << std::endl;
      }
    }
  }
  std::cout << "MIP: " << ncl << " clusters, last at y = "
            << yLast * 1.e4 << " um; " << nPrimary
            << " primary e-h pairs, " << nTotal
            << " electrons after multiplication -> counting gain "
            << double(nTotal) / double(nPrimary) << std::endl;
  std::cout << "[timer] MIP gain-ON pass: " << ElapsedS(tMipOnStart)
            << " s (" << primaries.size() << " primaries)" << std::endl;

  double qInt = 0., iMax = 0.;
  for (unsigned int i = 0; i < 800; ++i) {
    const double s = sensor.GetSignal("pad", i);
    qInt += s;
    if (std::abs(s) > std::abs(iMax)) iMax = s;
  }
  std::cout << "signal check (gain ON): peak bin = " << iMax
            << ", integral (arb.) = " << qInt << std::endl;

  // no-gain reference: same primaries, multiplication disabled
  sensor.ClearSignal();
  aval.EnableMultiplication(false);
  const auto tMipOffStart = Clock::now();
  for (const auto& p : primaries) {
    aval.AvalancheElectronHole(p[0], p[1], p[2], p[3]);
  }
  std::cout << "[timer] MIP gain-OFF pass: " << ElapsedS(tMipOffStart)
            << " s" << std::endl;
  double qIntNoGain = 0., iMaxNoGain = 0.;
  for (unsigned int i = 0; i < 800; ++i) {
    const double s = sensor.GetSignal("pad", i);
    qIntNoGain += s;
    if (std::abs(s) > std::abs(iMaxNoGain)) iMaxNoGain = s;
  }
  std::cout << "signal check (gain OFF, reference): peak bin = "
            << iMaxNoGain << ", integral (arb.) = " << qIntNoGain
            << std::endl;
  const double chargeGain = std::abs(qIntNoGain) > 0.
      ? qInt / qIntNoGain : 0.;
  std::cout << "charge-based gain (integral with gain / integral without) "
            << "= " << chargeGain << std::endl;
  aval.EnableMultiplication(true);

  // restore gain-ON signal for the plots below
  sensor.ClearSignal();
  const auto tMipRestoreStart = Clock::now();
  for (const auto& p : primaries) {
    aval.AvalancheElectronHole(p[0], p[1], p[2], p[3]);
  }
  std::cout << "[timer] MIP restore pass: " << ElapsedS(tMipRestoreStart)
            << " s" << std::endl;
  std::cout << "[stepfn aval] calls=" << avalFnCalls.load()
            << " fine=" << avalFnFine.load() << " coarse="
            << avalFnCoarse.load() << std::endl;
  std::cout << "[timer] total run so far: " << ElapsedS(tRunStart)
            << " s" << std::endl;

  driftView.Plot2d();

  TCanvas c("c", "", 800, 600);
  ViewSignal vs;
  vs.SetSensor(&sensor);
  vs.SetCanvas(&c);
  vs.PlotSignal("pad");
  c.SaveAs("signal_pad.pdf");

  double t0 = 0.;
  double dt = 0.05;
  TCanvas canvas("c", "", 1400, 600);
  canvas.Divide(2, 1);
  auto pad1 = canvas.cd(1);
  auto pad2 = canvas.cd(2);
  driftView.SetCanvas((TPad*)canvas.cd(1));
  vs.SetCanvas((TPad*)canvas.cd(2));

  const std::size_t nFrames = 120;
  const auto tMovieStart = Clock::now();
  for (std::size_t i = 0; i < 120; ++i) {
    driftView.Clear();
    aval.SetTimeWindow(t0, t0 + dt);
    aval.ResumeAvalanche();
    driftView.Plot2d(true, true);
    vs.PlotSignal("pad");
    gSystem->ProcessEvents();
    constexpr bool gif = true;
    if (!gif) {
      char filename[50];
      snprintf(filename, 50, "frames/frame_%03zu.png", i);
      canvas.SaveAs(filename);
    } else {
      if (i == nFrames - 1) {
        canvas.Print("planar_movie.gif++");
      } else {
        canvas.Print("planar_movie.gif+3");
      }
    }
    t0 += dt;
  }
  std::cout << "[timer] movie loop (120 frames): " << ElapsedS(tMovieStart)
            << " s" << std::endl;
  std::cout << "[timer] TOTAL RUNTIME: " << ElapsedS(tRunStart)
            << " s" << std::endl;
  std::cout << "Done.\n";
  
  // TCanvas b("b", "", 800, 600);
  // wfieldView.PlotContourWeightingField("strip", "v");
  // b.SaveAs("weighting_field.pdf");

  app.Run();
}