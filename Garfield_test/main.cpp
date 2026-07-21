#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <vector>

#include "Garfield/ComponentTcad2d.hh"
#include "Garfield/ComponentConstant.hh"
#include "Garfield/GeometrySimple.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/SolidBox.hh"
#include "Garfield/TrackHeed.hh"

#include "lgad_tools.hh"

using namespace Garfield;

struct RunConfig {
  std::string file = "/home/ahaines561/HEP/MAS/Silvaco_dat/"
                        // "diode.sta";
                        // "lgad150V.sta";
                        // "lgad180V.sta";
                        // "lgad190V.sta";
                        "aclgad.sta";
  std::string outDir = "output_files_aclgad/";

  double biasVOverride = std::numeric_limits<double>::quiet_NaN();
  double xTrackUm = 250.;  // MIP track position (250 = under strip1)
  bool doWeightingDump = true; // DumpWeightingField
  bool doConvergenceScan = false;  // RunConvergenceScan
  bool doModelComparison = false; // RunModelComparison
  bool doMIP = true;
  std::string model = "okuto";
  double driftWindowNs = 6.; // must exceed the 4ns signal window
  double fineStepNm = 20.; // step size inside the gain-layer band
  double bulkStepNm = 250.;  // step size everywhere else
  double fineBandHalfWidthUm = 2.5;  // band = [yGain-this, yGain+this]
};

int main() {
  RunConfig cfg;
  const auto tRunStart = Clock::now();

  double biasV = cfg.biasVOverride;
  if (std::isnan(biasV)) biasV = ParseBiasFromFilename(cfg.file);
  const std::string biasLabel = FormatBias(biasV);

  std::filesystem::create_directories(cfg.outDir);
  std::cout << "bias = " << biasLabel << (std::isnan(biasV) ? "" : " V")
            << ", output directory = " << cfg.outDir << std::endl;

  ComponentTcad2d cmp;
  if (!cmp.InitialiseSilvaco(cfg.file)) return 1;

  MediumSilicon si;
  si.SetTemperature(300);
  if (cfg.model == "massey") {
    si.SetImpactIonisationModelMassey();
  } else if (cfg.model == "grant") {
    si.SetImpactIonisationModelGrant();
  } else if (cfg.model == "okuto") {
    si.SetImpactIonisationModelOkutoCrowell();
  } else {
    si.SetImpactIonisationModelVanOverstraetenDeMan();
  }
  std::cout << "impact-ionisation model: "
            << (cfg.model == "massey" || cfg.model == "grant" ||
                        cfg.model == "okuto"
                    ? cfg.model : "vodm (van Overstraeten-de Man)")
            << std::endl;
  for (unsigned int i = 0; i < cmp.GetNumberOfRegions(); ++i) {
    cmp.SetMedium(i, &si);
  }
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
  ScanValidity(cmp, bx0, bx1, by0, by1, 400, 400,
              cfg.outDir + "validity_map.csv");

  double x0 = cfg.xTrackUm * 1.e-4 + 0.13e-4;
  if (x0 <= bx0 || x0 >= bx1) {
    std::cout << "requested x = " << cfg.xTrackUm
              << " um is outside the map; using the device centre."
              << std::endl;
    x0 = 0.5 * (bx0 + bx1) + 0.13e-4;
  }
  std::cout << "probe/track line at x = " << x0 * 1.e4 << " um"
            << std::endl;

  const auto prof = ScanFieldProfile(cmp, x0, by0, by1, 400,
                                     cfg.outDir + (biasLabel == "NA"
                                         ? "profile.csv"
                                         : "profile_" + biasLabel + "V.csv"));
  if (!prof.valid) return 1;
  const double yTop = prof.yTop;
  const double yBot = prof.yBot;
  const double d = prof.d;
  const double yGain = prof.yGain;
  const double eMax = prof.eMax;

  // 2D search for a stronger field than the 1D probe line found. Seeded
  // with the probe line's own peak (not 0) so the line's result is a
  // valid answer if the grid search finds nothing better, and so the
  // "e > globalMaxField" check below is correct from the first iteration.
  double globalMaxField = eMax;
  double globalMaxX = x0;
  double globalMaxY = yGain;
  const auto tHotspotStart = Clock::now();
  for (int ix = 0; ix <= 60; ++ix) {
    const double x = bx0 + ((bx1 - bx0) * ix) / 60. + 0.077e-4;
    for (int iy = 0; iy <= 150; ++iy) {
      const double y = by0 + ((by1 - by0) * iy) / 150. + 0.0113e-4;
      double ex, ey, ez, v; int st; Medium* m = nullptr;
      cmp.ElectricField(x, y, 0., ex, ey, ez, v, m, st);
      if (st != 0) continue;
      const double e = std::sqrt(ex * ex + ey * ey);
      if (e > globalMaxField) { globalMaxField = e; globalMaxX = x; globalMaxY = y; }
    }
    for (int iy = 0; iy <= 120; ++iy) {
      const double y = yTop - 0.5e-4 + (3.5e-4 * iy) / 120.;
      double ex, ey, ez, v; int st; Medium* m = nullptr;
      cmp.ElectricField(x, y, 0., ex, ey, ez, v, m, st);
      if (st != 0) continue;
      const double e = std::sqrt(ex * ex + ey * ey);
      if (e > globalMaxField) { globalMaxField = e; globalMaxX = x; globalMaxY = y; }
    }
  }
  std::cout << "global max field " << globalMaxField << " V/cm at (x, y) = ("
            << globalMaxX * 1.e4 << ", " << globalMaxY * 1.e4 << ") um"
            << std::endl;
  std::cout << "[timer] hotspot scan: " << ElapsedS(tHotspotStart)
            << " s" << std::endl;

  // field dump over x = xTrack +- 5 um, full depth, for offline plotting
  {
    const auto tDumpStart = Clock::now();
    const double xLo = std::max(bx0, (cfg.xTrackUm - 5.) * 1.e-4);
    const double xHi = std::min(bx1, (cfg.xTrackUm + 5.) * 1.e-4);
    const double dxCm = 0.1e-4, dyCm = 0.05e-4;  // 0.1um cols, 0.05um rows
    const int nx = static_cast<int>(std::round((xHi - xLo) / dxCm));
    const int ny = static_cast<int>(std::round((by1 - by0) / dyCm));
    const auto r = DumpElectricField(cmp, xLo, xHi, by0, by1, nx, ny,
                                     cfg.outDir + "efield_window.csv");
    std::cout << "field dump: x = [" << xLo * 1.e4 << ", " << xHi * 1.e4
              << "] um -> efield_window.csv (" << r.nRows << " points, "
              << r.nBad << " invalid skipped)" << std::endl;
    std::cout << "  window max |E| = " << r.wMax << " V/cm at (x, y) = ("
              << r.wMaxX * 1.e4 << ", " << r.wMaxY * 1.e4 << ") um"
              << std::endl;
    std::cout << "[timer] field dump: " << ElapsedS(tDumpStart) << " s"
              << std::endl;
  }
  {
    const auto tFullStart = Clock::now();
    const auto r = DumpElectricField(cmp, bx0, bx1, by0, by1, 250, 250,
                                     cfg.outDir + "efield_full.csv");
    std::cout << "full-device dump -> efield_full.csv (" << r.nRows
              << " points, " << r.nBad << " invalid skipped)" << std::endl;
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
    std::ofstream fprim(cfg.outDir + (biasLabel == "NA" ? "primary_pairs.txt"
        : ("primary_pairs_" + biasLabel + "V.txt")));
    for (int it = 0; it < nHeedTracks; ++it) {
      htrack.NewTrack(x0, yTop + 0.03e-4, 0., 0., 0., 1., 0.);
      double xc, yc, zc, tc, ec, extra;
      int nc = 0;
      unsigned long nPrimary = 0;
      while (htrack.GetCluster(xc, yc, zc, tc, nc, ec, extra)) {
        nPrimary += nc;
      }
      fprim << nPrimary << "\n";
    }
    std::cout << "wrote " << nHeedTracks << " fast Heed tracks to "
              << "primary_pairs" << (biasLabel == "NA" ? "" : "_" +
                 biasLabel + "V") << ".txt" << std::endl;
    std::cout << "[timer] Heed primary-statistics: "
              << ElapsedS(tHeedStart) << " s" << std::endl;
  }

  // 3 strips from the mask: edge -> centre/half-width
  //   M1 30-70->50,20  M2 220-270->245,25  M3 430-470->450,20
  const std::vector<Strip> strips = {
    {"strip0",  50., 20.},
    {"strip1", 245., 25.},
    {"strip2", 450., 20.},
  };
  if (!StripsInsideMap(strips, bx0, bx1)) return 1;

  ComponentAnalyticField wcmp;
  wcmp.AddPlaneY(yTop, 1., "top");
  wcmp.AddPlaneY(yBot, 0., "back");
  for (const auto& s : strips) {
    const double xc = s.centerUm * 1.e-4, hw = s.halfWidthUm * 1.e-4;
    wcmp.AddStripOnPlaneY('z', yTop, xc - hw, xc + hw, s.label);
  }

  PrintWeightingSanity(wcmp, strips, yTop, yBot);

  Sensor sensor;
  sensor.AddComponent(&cmp);
  for (const auto& s : strips) sensor.AddElectrode(&wcmp, s.label);
  std::cout << "Sensor has " << sensor.GetNumberOfElectrodes()
            << " electrodes (expect " << strips.size() << ")" << std::endl;

  if (cfg.doWeightingDump) {
    DumpWeightingField(wcmp, strips, bx0, bx1, yTop, yBot,
                       cfg.outDir + "wfield_full.txt");
  }

  sensor.SetTimeWindow(0., 0.005, 800);
  sensor.SetArea(bx0, yTop + 0.02e-4, -5.e-4, bx1, yBot, 5.e-4);

  const double fineStepCm = cfg.fineStepNm * 1.e-7;   // nm -> cm
  const double bulkStepCm = cfg.bulkStepNm * 1.e-7;
  const double fineBandCm = cfg.fineBandHalfWidthUm * 1.e-4;
  const double yFineLo = yGain - fineBandCm, yFineHi = yGain + fineBandCm;

  double avalStepCm = fineStepCm;  // mutable -- captured by ref in the
                                    // step function set up below
  std::atomic<long long> avalCalls{0}, avalFine{0}, avalCoarse{0},
      avalPrinted{0};
  AvalancheMC aval;
  aval.SetSensor(&sensor);
  const std::size_t sizeCap = 5000;
  ConfigureAvalanche(aval, avalStepCm, bulkStepCm, yFineLo, yFineHi,
                     cfg.driftWindowNs, sizeCap, /*enableSignal=*/true,
                     /*multithreading=*/true, "aval", avalCalls, avalFine,
                     avalCoarse, avalPrinted);

  const double yInj = std::min(yGain + 5.e-4, 0.5 * (yGain + yBot));
  double ladderStepCm = fineStepCm;
  std::atomic<long long> ladderCalls{0}, ladderFine{0}, ladderCoarse{0},
      ladderPrinted{0};
  AvalancheMC avalLadder;
  avalLadder.SetSensor(&sensor);
  const std::size_t ladderCap = 5000;
  ConfigureAvalanche(avalLadder, ladderStepCm, bulkStepCm, yFineLo, yFineHi,
                     cfg.driftWindowNs, ladderCap, /*enableSignal=*/false,
                     /*multithreading=*/true, "ladder", ladderCalls,
                     ladderFine, ladderCoarse, ladderPrinted);

  if (cfg.doConvergenceScan) {
    RunConvergenceScan(avalLadder, ladderStepCm, x0, yInj, ladderCap,
                       cfg.outDir);
    std::cout << "[stepfn ladder] calls=" << ladderCalls.load()
              << " fine=" << ladderFine.load() << " coarse="
              << ladderCoarse.load() << std::endl;
  }

  if (cfg.doModelComparison) {
    RunModelComparison(avalLadder, si, x0, yInj, ladderStepCm, ladderCap,
                       cfg.outDir, biasLabel, eMax, globalMaxField);
    std::cout << "[stepfn ladder] calls=" << ladderCalls.load()
              << " fine=" << ladderFine.load() << " coarse="
              << ladderCoarse.load() << std::endl;
    // restore the configured MIP model -- RunModelComparison leaves si
    // set to whichever model it evaluated last
    if (cfg.model == "massey") {
      si.SetImpactIonisationModelMassey();
    } else if (cfg.model == "grant") {
      si.SetImpactIonisationModelGrant();
    } else if (cfg.model == "okuto") {
      si.SetImpactIonisationModelOkutoCrowell();
    } else {
      si.SetImpactIonisationModelVanOverstraetenDeMan();
    }
  }

  if (!cfg.doMIP) {
    std::cout << "doMIP=false -- nothing more to compute." << std::endl;
    return 0;
  }

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

  track.NewTrack(x0, yTop + 0.03e-4, 0., 0., 0., 1., 0.);
  double xc, yc, zc, tc, ec, extra;
  int nc = 0;
  unsigned long nPrimary = 0; // total e-h pairs, all clusters
  unsigned long ncl = 0;  // cluster count
  unsigned long nTotal = 0; // electrons after multiplication (gain-ON loop)
  unsigned long nPairsDone = 0; // gain-ON loop progress counter
  double yLast = 0.;

  std::vector<std::array<double, 4>> primaries;
  while (track.GetCluster(xc, yc, zc, tc, nc, ec, extra)) {
    ++ncl;
    yLast = yc;
    nPrimary += nc;
    for (int k = 0; k < nc; ++k) primaries.push_back({xc, yc, zc, tc});
  }
  std::cout << "MIP track: " << ncl << " clusters, last at y = "
            << yLast * 1.e4 << " um, " << nPrimary
            << " primary e-h pairs" << std::endl;

  sensor.ClearSignal();
  aval.EnableMultiplication(false);
  const auto tMipOffStart = Clock::now();
  {
    unsigned long nOffDone = 0;
    for (const auto& p : primaries) {
      aval.AvalancheElectronHole(p[0], p[1], p[2], p[3]);
      if (++nOffDone % 500 == 0) {
        std::cout << "  [OFF] pair " << nOffDone << "/" << primaries.size()
                  << "  [" << ElapsedS(tMipOffStart) << " s]" << std::endl;
      }
    }
  }
  std::cout << "[timer] MIP gain-OFF pass: " << ElapsedS(tMipOffStart)
            << " s" << std::endl;
  std::vector<std::vector<double>> sigOff(strips.size(),
                                          std::vector<double>(800, 0.));
  std::vector<double> qIntNoGain(strips.size(), 0.),
      iMaxNoGain(strips.size(), 0.);
  for (std::size_t k = 0; k < strips.size(); ++k) {
    for (unsigned int i = 0; i < 800; ++i) {
      const double s = sensor.GetSignal(strips[k].label, i);
      sigOff[k][i] = s;
      qIntNoGain[k] += s;
      if (std::abs(s) > std::abs(iMaxNoGain[k])) iMaxNoGain[k] = s;
    }
    std::cout << strips[k].label << " signal check (gain OFF): peak bin = "
              << iMaxNoGain[k] << ", integral (arb.) = " << qIntNoGain[k]
              << std::endl;
  }

  sensor.ClearSignal();
  aval.EnableMultiplication(true);
  std::ofstream fpairs(cfg.outDir + (biasLabel == "NA" ? "mip_pairs.csv"
      : ("mip_pairs_" + biasLabel + "V.csv")));
  fpairs << "x_um,y_um,ne\n";
  const auto tMipOnStart = Clock::now();
  for (const auto& p : primaries) {
    aval.AvalancheElectronHole(p[0], p[1], p[2], p[3]);
    std::size_t ne = 0, ni = 0;
    aval.GetAvalancheSize(ne, ni);
    nTotal += ne;
    fpairs << p[0] * 1.e4 << "," << p[1] * 1.e4 << "," << ne << "\n";
    if (++nPairsDone % 200 == 0) {
      std::cout << "  pair " << nPairsDone << "/" << primaries.size()
                << "  [" << ElapsedS(tMipOnStart) << " s elapsed]"
                << std::endl;
    }
  }
  fpairs.close();
  const double countingGain = double(nTotal) / double(nPrimary);
  std::cout << "MIP: " << nTotal << " electrons after multiplication -> "
            << "counting gain " << countingGain << std::endl;
  std::cout << "[timer] MIP gain-ON pass: " << ElapsedS(tMipOnStart)
            << " s (" << primaries.size() << " primaries)" << std::endl;

  std::vector<std::vector<double>> sigOn(strips.size(),
                                         std::vector<double>(800, 0.));
  std::vector<double> qInt(strips.size(), 0.), iMax(strips.size(), 0.);
  for (std::size_t k = 0; k < strips.size(); ++k) {
    for (unsigned int i = 0; i < 800; ++i) {
      const double s = sensor.GetSignal(strips[k].label, i);
      sigOn[k][i] = s;
      qInt[k] += s;
      if (std::abs(s) > std::abs(iMax[k])) iMax[k] = s;
    }
    const double stripGain = std::abs(qIntNoGain[k]) > 0.
        ? qInt[k] / qIntNoGain[k] : 0.;
    std::cout << strips[k].label << " signal check (gain ON): peak bin = "
              << iMax[k] << ", integral (arb.) = " << qInt[k]
              << ", charge-based gain = " << stripGain << std::endl;
  }

  {
    double qTotalOn = 0.;
    for (const auto& q : qInt) qTotalOn += std::abs(q);
    std::cout << "charge sharing @ x0=" << x0 * 1.e4 << " um: ";
    for (std::size_t k = 0; k < strips.size(); ++k) {
      const double frac = qTotalOn > 0. ? std::abs(qInt[k]) / qTotalOn : 0.;
      std::cout << strips[k].label << "=" << frac << " ";
    }
    std::cout << std::endl;

    std::ifstream fcstest(cfg.outDir + "charge_sharing.csv");
    const bool csExists = fcstest.good();
    fcstest.close();
    std::ofstream fcs(cfg.outDir + "charge_sharing.csv", std::ios::app);
    if (!csExists) {
      fcs << "bias,model,x0_um";
      for (const auto& s : strips) fcs << "," << s.label << "_frac";
      for (const auto& s : strips) fcs << "," << s.label << "_qOn_arb";
      for (const auto& s : strips) fcs << "," << s.label << "_qOff_arb";
      fcs << "\n";
    }
    fcs << biasLabel << "," << cfg.model << "," << x0 * 1.e4;
    for (std::size_t k = 0; k < strips.size(); ++k) {
      const double frac = qTotalOn > 0. ? std::abs(qInt[k]) / qTotalOn : 0.;
      fcs << "," << frac;
    }
    for (std::size_t k = 0; k < strips.size(); ++k) fcs << "," << qInt[k];
    for (std::size_t k = 0; k < strips.size(); ++k)
      fcs << "," << qIntNoGain[k];
    fcs << "\n";
    std::cout << "appended charge-sharing row to charge_sharing.csv"
              << std::endl;
  }

  double qIntSumOn = 0., qIntSumOff = 0.;
  for (std::size_t k = 0; k < strips.size(); ++k) {
    qIntSumOn += qInt[k];
    qIntSumOff += qIntNoGain[k];
  }
  const double chargeGain = std::abs(qIntSumOff) > 0.
      ? qIntSumOn / qIntSumOff : 0.;
  std::cout << "device charge-based gain (sum over all strips, ON/OFF) = "
            << chargeGain << std::endl;
  std::size_t kPeak = 0;
  for (std::size_t k = 1; k < strips.size(); ++k) {
    if (std::abs(iMax[k]) > std::abs(iMax[kPeak])) kPeak = k;
  }

  {
    std::ofstream fsig(cfg.outDir + (biasLabel == "NA" ? "signal_pad.csv"
        : ("signal_" + biasLabel + "V.csv")));
    fsig << "bin,t_ns";
    for (const auto& s : strips) fsig << "," << s.label << "_ON," << s.label
                                       << "_OFF";
    fsig << "\n";
    for (unsigned int i = 0; i < 800; ++i) {
      fsig << i << "," << (i + 0.5) * 0.005;
      for (std::size_t k = 0; k < strips.size(); ++k) {
        fsig << "," << sigOn[k][i] << "," << sigOff[k][i];
      }
      fsig << "\n";
    }
  }
  {
    std::ifstream ftest(cfg.outDir + "mip_summary.csv");
    const bool exists = ftest.good();
    ftest.close();
    std::ofstream fsum(cfg.outDir + "mip_summary.csv", std::ios::app);
    if (!exists) {
      fsum << "bias,model,nClusters,nPrimary,countingGain,chargeGain,"
           << "peakOn,peakOff,intOn,intOff,elapsed_s\n";
    }
    fsum << biasLabel << "," << cfg.model << "," << ncl << "," << nPrimary
         << "," << countingGain << "," << chargeGain << ","
         << iMax[kPeak] << "," << iMaxNoGain[kPeak] << "," << qIntSumOn
         << "," << qIntSumOff << "," << ElapsedS(tMipOnStart) << "\n";
    std::cout << "appended MIP summary row (bias=" << biasLabel
              << ", model=" << cfg.model << ", peak strip="
              << strips[kPeak].label << ") to mip_summary.csv" << std::endl;
  }
  std::cout << "[stepfn aval] calls=" << avalCalls.load()
            << " fine=" << avalFine.load() << " coarse="
            << avalCoarse.load() << std::endl;
  std::cout << "[timer] total run so far: " << ElapsedS(tRunStart)
            << " s" << std::endl;

  std::cout << "[timer] TOTAL RUNTIME: " << ElapsedS(tRunStart)
            << " s" << std::endl;
  std::cout << "Done.\n";
}