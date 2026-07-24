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

int main(int argc, char* argv[]) {
  const auto tRunStart = Clock::now();

  // CLI Arguments
  const std::string file = argc > 1 ? argv[1] : "/home/ahaines561/HEP/MAS/Silvaco_dat/lgad190V.sta";
  const int nMips = argc > 2 ? std::atoi(argv[2]) : 200; //num of MIPs
  const std::string biasV = "180V";
  const std::string outDir = "output_files_1000MIPs/";
  std::filesystem::create_directories(outDir);

  std::cout << "Starting overlay run with " << nMips << " MIPs." << std::endl;
  std::cout << "Loading field map: " << file << std::endl;

  ComponentTcad2d cmp;
  if (!cmp.InitialiseSilvaco(file)) return 1;

  MediumSilicon si;
  si.SetTemperature(293.15);
  si.SetImpactIonisationModelOkutoCrowell();
  cmp.SetMedium(0, &si);
  cmp.SetRangeZ(-5.e-4, 5.e-4);

  // Get Bounding Box
  double bx0 = 0., by0 = 0., bz0 = 0., bx1 = 0., by1 = 0., bz1 = 0.;
  if (!cmp.GetBoundingBox(bx0, by0, bz0, bx1, by1, bz1)) {
    std::cerr << "Could not get bounding box.\n";
    return 1;
  }

  double yTop = 1., yBot = -1., eMax = 0., yGain = 0.;
  const double kActiveFieldMinVcm = 100.;
  // const double x0 = 250.13e-4;
    const double x0 = 20.13e-4;
  for (int i = 0; i <= 400; ++i) {
    const double y = by0 + ((by1 - by0) * i) / 400.;
    double ex, ey, ez, v; int st; Medium* m = nullptr;
    cmp.ElectricField(x0, y, 0., ex, ey, ez, v, m, st);
    if (st != 0) continue;
    const double e = std::sqrt(ex * ex + ey * ey);
    if (e < kActiveFieldMinVcm) continue;
    if (yTop > yBot) yTop = y;
    yBot = y;
    if (e > eMax) { eMax = e; yGain = y; }
  }
  const double d = yBot - yTop;
  std::cout << "Track at x = " << x0 * 1.e4 << " um" << std::endl;
  std::cout << "Active silicon thickness: " << d * 1.e4 << " um" << std::endl;
  std::cout << "Peak field " << eMax << " V/cm at y = " << yGain * 1.e4 << " um" << std::endl;

  // strips
  // const std::vector<Strip> strips = {
  //   {"strip0",  50., 20.},
  //   {"strip1", 245., 25.},
  //   {"strip2", 450., 20.},
  // };
  const std::vector<Strip> strips = {
    {"anode",  22.5, 22.5},
    {"cathode", 77.5, 22.5},
  };
  if (!StripsInsideMap(strips, bx0, bx1)) return 1;
  const std::size_t nStrips = strips.size();

  ComponentAnalyticField wcmp;
  wcmp.AddPlaneY(yTop, 1., "top");
  wcmp.AddPlaneY(yBot, 0., "back");
  for (const auto& s : strips) {
    const double xc = s.centerUm * 1.e-4, hw = s.halfWidthUm * 1.e-4;
    wcmp.AddStripOnPlaneY('z', yTop, xc - hw, xc + hw, s.label);
  }

  PrintWeightingSanity(wcmp, strips, yTop, yBot);

  const double yFine = yGain + 2.5e-4;

  // HEED
  std::cout << "\nPhase 1: Generating " << nMips << " tracks sequentially (HEED)..." << std::endl;

  SolidBox box(0.5 * (bx0 + bx1), 0.5 * (yTop + yBot), 0., 0.5 * (bx1 - bx0), 0.5 * d, 5.e-4);
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

  // Store primary coordinates
  std::vector<std::vector<std::array<double, 4>>> allPrimaries(nMips);

  for (int iTrk = 0; iTrk < nMips; ++iTrk) {
    track.NewTrack(x0, yTop + 0.03e-4, 0., 0., 0., 1., 0.);
    double xc, yc, zc, tc, ec, extra;
    int nc = 0;
    while (track.GetCluster(xc, yc, zc, tc, nc, ec, extra)) {
      for (int k = 0; k < nc; ++k) allPrimaries[iTrk].push_back({xc, yc, zc, tc});
    }
  }
  std::cout << "Phase 1 complete. Proceeding to Avalanche simulation." << std::endl;

  // process the biggest (slowest) tracks first: with dynamic scheduling,
  // a large track that lands late on an already-busy thread becomes a
  // lone straggler holding up the whole batch at the end. This only
  // reorders which track index a thread works on, not the physics.
  std::vector<int> trackOrder(nMips);
  for (int i = 0; i < nMips; ++i) trackOrder[i] = i;
  std::sort(trackOrder.begin(), trackOrder.end(), [&](int a, int b) {
    return allPrimaries[a].size() > allPrimaries[b].size();
  });

  // OpenMP avalanche loop -- [strip][track][bin]
  std::vector<std::vector<std::vector<double>>> allSignals(
      nStrips, std::vector<std::vector<double>>(nMips,
                                                std::vector<double>(800, 0.)));
  std::vector<std::vector<std::vector<double>>> allSignalsOff(
      nStrips, std::vector<std::vector<double>>(nMips,
                                                std::vector<double>(800, 0.)));

  std::cout << "\nPhase 2: Starting " << nMips << " MIP avalanches across 8 cores..." << std::endl;
  std::cout << "(Silencing Garfield++ OpenMP oversubscription warnings...)" << std::endl;

  std::ofstream devNull("/dev/null");
  std::streambuf* oldCerr = std::cerr.rdbuf(devNull.rdbuf());

  // gain-OFF reference pass: same primaries, no multiplication. Far
  // cheaper than the ON pass (no cascade), so it runs first.
  std::cout << "\nPhase 2a: gain-OFF reference pass..." << std::endl;
  const auto tOffStart = Clock::now();

  #pragma omp parallel
  {
    Sensor localSensor;
    localSensor.AddComponent(&cmp);
    for (const auto& s : strips) localSensor.AddElectrode(&wcmp, s.label);
    localSensor.SetTimeWindow(0., 0.005, 800);
    localSensor.SetArea(bx0, yTop + 0.02e-4, -5.e-4, bx1, yBot, 5.e-4);

    AvalancheMC localAval;
    localAval.EnableMultithreading(false);
    localAval.SetSensor(&localSensor);
    localAval.EnableSignalCalculation();
    localAval.EnableMultiplication(false);
    localAval.SetTimeWindow(0., 6.);
    localAval.EnableAvalancheSizeLimit(5000);
    localAval.SetStepDistanceFunction([yFine](double x, double y, double z) {
      if (y < yFine) return 5.e-6 ;
      return 2.5e-5;
    });

    #pragma omp for schedule(dynamic, 1)
      for (int k = 0; k < nMips; ++k) {
        const int iTrk = trackOrder[k];
        localSensor.ClearSignal();
        for (const auto& p : allPrimaries[iTrk]) {
          localAval.AvalancheElectronHole(p[0], p[1], p[2], p[3]);
        }
        for (std::size_t is = 0; is < nStrips; ++is) {
          for (unsigned int i = 0; i < 800; ++i) {
            allSignalsOff[is][iTrk][i] =
                localSensor.GetSignal(strips[is].label, i);
          }
        }
        #pragma omp critical
        {
          std::cout << "  [OFF] Track " << iTrk + 1 << "/" << nMips
                    << " (pairs=" << allPrimaries[iTrk].size() << ")"
                    << " | Elapsed: " << ElapsedS(tOffStart) << " s"
                    << std::endl;
        }
      }
  }
  std::cout << "Phase 2a complete: " << ElapsedS(tOffStart) << " s"
            << std::endl;

  // gain-ON pass
  std::cout << "\nPhase 2b: gain-ON pass..." << std::endl;
  const auto tOverlayStart = Clock::now();

  #pragma omp parallel
  {
    // Thread-local sensor and avalanche objects
    Sensor localSensor;
    localSensor.AddComponent(&cmp);
    for (const auto& s : strips) localSensor.AddElectrode(&wcmp, s.label);
    localSensor.SetTimeWindow(0., 0.005, 800);
    localSensor.SetArea(bx0, yTop + 0.02e-4, -5.e-4, bx1, yBot, 5.e-4);

    AvalancheMC localAval;
    localAval.EnableMultithreading(false); // Outer loop controls threads now
    localAval.SetSensor(&localSensor);
    localAval.EnableSignalCalculation();
    localAval.EnableMultiplication(true);
    // signal window is 4ns (800 x 0.005); without a drift bound a carrier
    // in a near-zero-field pocket (e.g. the entry layer) drifts forever
    localAval.SetTimeWindow(0., 6.);
    localAval.EnableAvalancheSizeLimit(5000);
    localAval.SetStepDistanceFunction([yFine](double x, double y, double z) {
      if (y < yFine) return 5.e-6 ;
      return 2.5e-5;
    });

    #pragma omp for schedule(dynamic, 1)
      for (int k = 0; k < nMips; ++k) {
        const int iTrk = trackOrder[k];
        localSensor.ClearSignal();

        int nOk = 0, nFail = 0;
        for (const auto& p : allPrimaries[iTrk]) {
          bool ok = localAval.AvalancheElectronHole(p[0], p[1], p[2], p[3]);
          if (ok) ++nOk; else ++nFail;
        }

        for (std::size_t is = 0; is < nStrips; ++is) {
          for (unsigned int i = 0; i < 800; ++i) {
            allSignals[is][iTrk][i] =
                localSensor.GetSignal(strips[is].label, i);
          }
        }

        #pragma omp critical
        {
          std::cout << "  [ON] Track " << iTrk + 1 << "/" << nMips
                    << " pairs=" << allPrimaries[iTrk].size()
                    << " | Elapsed: " << ElapsedS(tOverlayStart) << " s"
                    << std::endl;
        }
      }
  }

  std::cerr.rdbuf(oldCerr);

  // one CSV per strip per pass, same format as before so the existing
  // notebook loader works unchanged
  for (std::size_t is = 0; is < nStrips; ++is) {
    for (int pass = 0; pass < 2; ++pass) {
      const auto& sig = (pass == 0) ? allSignals[is] : allSignalsOff[is];
      const std::string tag = (pass == 0) ? "" : "_noGain";
      const std::string path =
          outDir + "signal_overlay" + tag + "_" + strips[is].label + biasV
          + ".csv";
      std::ofstream f(path);
      f << "time_ns";
      for (int iTrk = 0; iTrk < nMips; ++iTrk) f << ",trk" << iTrk;
      f << "\n";
      f << "nPairs";
      for (int iTrk = 0; iTrk < nMips; ++iTrk) {
        f << "," << allPrimaries[iTrk].size();
      }
      f << "\n";
      for (unsigned int i = 0; i < 800; ++i) {
        f << (i + 0.5) * 0.005;
        for (int iTrk = 0; iTrk < nMips; ++iTrk) f << "," << sig[iTrk][i];
        f << "\n";
      }
      std::cout << "wrote " << path << std::endl;
    }
  }

  std::cout << "TOTAL RUNTIME: " << ElapsedS(tRunStart) << " s" << std::endl;
  return 0;
}