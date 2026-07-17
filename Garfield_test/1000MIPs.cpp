#include <cmath>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <array>
#include <chrono>
#include <iostream>
#include <filesystem>
#include <iomanip>
#include <omp.h>

#include "Garfield/ComponentTcad2d.hh"
#include "Garfield/ComponentConstant.hh"
#include "Garfield/MediumSilicon.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/GeometrySimple.hh"
#include "Garfield/SolidBox.hh"
#include "Garfield/TrackHeed.hh"
#include "Garfield/AvalancheMC.hh"

using namespace Garfield;
using Clock = std::chrono::steady_clock;

double ElapsedS(const Clock::time_point& t0) {
  return std::chrono::duration<double>(Clock::now() - t0).count();
}

int main(int argc, char* argv[]) {
  const auto tRunStart = Clock::now();

  // CLI Arguments
  const std::string file = argc > 1 ? argv[1] : "/home/ahaines561/HEP/MAS/Silvaco_dat/lgad150V.sta";
  const int nMips = argc > 2 ? std::atoi(argv[2]) : 200; //num of MIPs
  const std::string biasV = "150V";
  const std::string outDir = "output_files/";
  std::filesystem::create_directories(outDir);

  std::cout << "Starting overlay run with " << nMips << " MIPs." << std::endl;
  std::cout << "Loading field map: " << file << std::endl;

  ComponentTcad2d cmp;
  if (!cmp.InitialiseSilvaco(file)) return 1;

  MediumSilicon si;
  si.SetTemperature(300);
  si.SetImpactIonisationModelOkutoCrowell();

  cmp.SetMedium("3", &si);
  cmp.SetRangeZ(-5.e-4, 5.e-4);

  // Get Bounding Box
  double bx0 = 0., by0 = 0., bz0 = 0., bx1 = 0., by1 = 0., bz1 = 0.;
  if (!cmp.GetBoundingBox(bx0, by0, bz0, bx1, by1, bz1)) {
    std::cerr << "Could not get bounding box.\n";
    return 1;
  }

  // Find active silicon depth (yTop to yBot)
  double yTop = 1., yBot = -1., eMax = 0., yGain = 0.;
  const double x0 = 20.13e-4; // Central track position
  for (int i = 0; i <= 400; ++i) {
    const double y = by0 + ((by1 - by0) * i) / 400.;
    double ex, ey, ez, v; int st; Medium* m = nullptr;
    cmp.ElectricField(x0, y, 0., ex, ey, ez, v, m, st);
    if (st != 0) continue; 
    if (yTop > yBot) yTop = y;
    yBot = y;
    const double e = std::sqrt(ex * ex + ey * ey);
    if (e > eMax) { eMax = e; yGain = y; }
  }
  const double d = yBot - yTop;
  std::cout << "Active silicon thickness: " << d * 1.e4 << " um" << std::endl;
  std::cout << "Peak field " << eMax << " V/cm at y = " << yGain * 1.e4 << " um" << std::endl;

  // Linear Weighting Field Approximation
  ComponentConstant wcmp;
  wcmp.SetArea(bx0, yTop, -5.e-4, bx1, yBot, 5.e-4);
  wcmp.SetMedium(&si);
  wcmp.SetElectricField(0., 0., 0.);
  wcmp.SetWeightingField(0., 1. / d, 0., "pad");
  wcmp.SetWeightingPotential(0.5 * (bx0 + bx1), yTop, 0., 1.);

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

  // OpenMP avalanche loop
  std::vector<std::vector<double>> allSignals(nMips, std::vector<double>(800, 0.));
  std::cout << "\nPhase 2: Starting " << nMips << " MIP avalanches across 8 cores..." << std::endl;
  std::cout << "(Silencing Garfield++ OpenMP oversubscription warnings...)" << std::endl;
  const auto tOverlayStart = Clock::now();

  std::ofstream devNull("/dev/null");
  std::streambuf* oldCerr = std::cerr.rdbuf(devNull.rdbuf());

  #pragma omp parallel
  {
    // Thread-local sensor and avalanche objects
    Sensor localSensor;
    localSensor.AddComponent(&cmp);
    localSensor.AddElectrode(&wcmp, "pad");
    localSensor.SetTimeWindow(0., 0.005, 800);
    localSensor.SetArea(bx0, yTop + 0.02e-4, -5.e-4, bx1, yBot, 5.e-4);

    AvalancheMC localAval;
    localAval.EnableMultithreading(false); // Outer loop controls threads now
    localAval.SetSensor(&localSensor);
    localAval.EnableSignalCalculation();
    localAval.EnableMultiplication(true);
    localAval.SetStepDistanceFunction([yFine](double x, double y, double z) {
      if (y < yFine) return 2.e-6; 
      return 2.5e-5;
    });

    #pragma omp for schedule(dynamic, 1)
      for (int iTrk = 0; iTrk < nMips; ++iTrk) {
        localSensor.ClearSignal();

        int nOk = 0, nFail = 0;
        for (const auto& p : allPrimaries[iTrk]) {
          bool ok = localAval.AvalancheElectronHole(p[0], p[1], p[2], p[3]);
          if (ok) ++nOk; else ++nFail;
        }

        double sigSum = 0.;
        for (unsigned int i = 0; i < 800; ++i) {
          allSignals[iTrk][i] = localSensor.GetSignal("pad", i);
          sigSum += std::abs(allSignals[iTrk][i]);
        }

        #pragma omp critical
        {
          std::cout << "  Track " << iTrk + 1 << "/" << nMips
                    << " pairs=" << allPrimaries[iTrk].size()
                    // << " ok=" << nOk << " fail=" << nFail
                    // << " endpointsE=" << localAval.GetNumberOfElectronEndpoints()
                    // << " endpointsH=" << localAval.GetHoles().size()
                    // << " |signal|_sum=" << sigSum
                    // << " | Thread " << omp_get_thread_num()
                    << " | Elapsed: " << ElapsedS(tOverlayStart) << " s"
                    << std::endl;
        }
      }
  }

  std::cerr.rdbuf(oldCerr);

  //csv
  std::ofstream fsig(outDir + "signal_overlay" + biasV + ".csv");
  fsig << "time_ns";
  for (int iTrk = 0; iTrk < nMips; ++iTrk) fsig << ",trk" << iTrk;
  fsig << "\n";

  for (unsigned int i = 0; i < 800; ++i) {
    fsig << (i + 0.5) * 0.005;
    for (int iTrk = 0; iTrk < nMips; ++iTrk) {
      fsig << "," << allSignals[iTrk][i];
    }
    fsig << "\n";
  }
  
  std::cout << "\nOverlay CSV written to " << outDir << "signal_overlay" + biasV + ".csv" << std::endl;
  std::cout << "TOTAL RUNTIME: " << ElapsedS(tRunStart) << " s" << std::endl;
  return 0;
}