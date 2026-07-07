#include <cmath>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <iostream>
#include <TApplication.h>
#include <TCanvas.h>

#include "Garfield/ComponentTcad2d.hh"
#include "Garfield/ComponentConstant.hh"
#include "Garfield/MediumSilicon.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/GeometrySimple.hh"
#include "Garfield/SolidBox.hh"
#include "Garfield/TrackHeed.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/ViewSignal.hh"

using namespace Garfield;

int main(int argc, char* argv[]) {
  const bool doMIP = argc > 2 ? std::atoi(argv[2]) != 0 : true;
  const std::string model = argc > 3 ? argv[3] : "vodm";
  const std::string file = argc > 1 ? argv[1]
      : "/home/ahaines561/HEP/MAS/Silvaco_dat/lgad.sta";
  int rootArgc = 1;
  char* rootArgv[] = {argv[0], nullptr};
  TApplication app("app", &rootArgc, rootArgv);


  ComponentTcad2d cmp;
  if (!cmp.InitialiseSilvaco(file)) return 1;

  MediumSilicon si;
  si.SetTemperature(293.15);
  // Impact-ionisation model, selected on the command line
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

  // device geometry
  double bx0 = 0., by0 = 0., bz0 = 0., bx1 = 0., by1 = 0., bz1 = 0.;
  if (!cmp.GetBoundingBox(bx0, by0, bz0, bx1, by1, bz1)) {
    std::cerr << "Could not get the bounding box from the field map;\n"
              << "set the device extent manually in the macro.\n";
    return 1;
  }
  std::cout << "map extent: x = [" << bx0 * 1.e4 << ", " << bx1 * 1.e4
            << "] um, y = [" << by0 * 1.e4 << ", " << by1 * 1.e4
            << "] um" << std::endl;
  const double x0 = 0.5 * (bx0 + bx1) + 0.13e-4;

  // Scans the full bounding box
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
  const double d = yBot - yTop;   // drift thickness
  std::cout << "active silicon: y = [" << yTop * 1.e4 << ", "
            << yBot * 1.e4 << "] um (d = " << d * 1.e4 << " um)"
            << std::endl;
  std::cout << "Peak field " << eMax << " V/cm at y = " << yGain * 1.e4
            << " um  (gain layer)" << std::endl;
  if (eMax < 2.5e5) {
    std::cout << "NOTE: peak field < 250 kV/cm -- expect gain near 1 "
              << "at this bias." << std::endl;
  }

  // Weighting field pproximation
  ComponentConstant wcmp;
  wcmp.SetArea(bx0, yTop, -5.e-4, bx1, yBot, 5.e-4);
  wcmp.SetMedium(&si);
  wcmp.SetElectricField(0., 0., 0.);
  wcmp.SetWeightingField(0., 1. / d, 0., "pad");
  wcmp.SetWeightingPotential(0., 0., 0., 1.);
  {
    double wx = 0., wy = 0., wz = 0.;
    wcmp.WeightingField(125.e-4, 25.e-4, 0., wx, wy, wz, "pad");
    std::cout << "Ew at centre = (" << wx << ", " << wy << ", " << wz
              << ")  [expect (0, 200, 0)]" << std::endl;
    std::cout << "wpot: top = "
              << wcmp.WeightingPotential(125.e-4, 0., 0., "pad")
              << ", mid = "
              << wcmp.WeightingPotential(125.e-4, 25.e-4, 0., "pad")
              << ", back = "
              << wcmp.WeightingPotential(125.e-4, 50.e-4, 0., "pad")
              << "  [expect 1, 0.5, 0]" << std::endl;
  }

  Sensor sensor;
  sensor.AddComponent(&cmp);
  sensor.AddElectrode(&wcmp, "pad");
  sensor.SetTimeWindow(0., 0.005, 800);
  sensor.SetArea(bx0, yTop + 0.02e-4, -5.e-4, bx1, yBot - 0.5e-4, 5.e-4);

  AvalancheMC aval;
  aval.SetSensor(&sensor);
  aval.EnableSignalCalculation();
  // the CMOS-LGAD calibration (arXiv:2505.05974) needed
  // 0.1 ps steps for TCAD-level gain agreement.
  aval.SetTimeSteps(2.e-3); // 2 ps
  // use 5.e-4 or finer for quotable gain number
  const std::size_t sizeCap = 20000;
  const double yInj = std::min(yGain + 5.e-4, 0.5 * (yGain + yBot));
  AvalancheMC avalLadder;
  avalLadder.SetSensor(&sensor);
  avalLadder.EnableAvalancheSizeLimit(sizeCap);
  aval.EnableAvalancheSizeLimit(sizeCap);
  const double steps[] = {2.e-3};
  std::vector<std::size_t> sizes;
  std::size_t nCapped = 0;
  std::cout << "# dt[ps]   G_e (mean+-sem, N)      G_eh (mean+-sem, N)\n";
  for (const double dt : steps) {
    avalLadder.SetTimeSteps(dt);
    double res[2][2] = {{0., 0.}, {0., 0.}};
    int nDone[2] = {0, 0};
    const int nWant[2] = {200, 300};
    for (int mode = 0; mode < 2; ++mode) {
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
        if (mode == 1) {
          sizes.push_back(ne);
          if (ne >= sizeCap) ++nCapped;
        }
      }
      // Print each mode's result as soon as it completes
      // wedged avalanche in the other mode cannot hide it
      const int n = nDone[mode];
      const double mean = res[mode][0] / n;
      const double var = res[mode][1] / n - mean * mean;
      const double sem = std::sqrt(var > 0. ? var / n : 0.);
      std::cout << "dt = " << dt * 1.e3 << " ps   "
                << (mode == 0 ? "G_e  = " : "G_eh = ")
                << mean << " +- " << sem << " (N=" << n << ")"
                << std::endl;
    }
  }

  if (!sizes.empty()) {
    std::size_t mx = 0, nBig = 0;
    double mean = 0.;
    for (const auto s : sizes) mean += s;
    mean /= sizes.size();
    for (const auto s : sizes) {
      if (s > mx) mx = s;
      if (s > 5. * mean) ++nBig;
    }
    std::cout << "e+h avalanche tail: max = " << mx << ", "
              << nBig << "/" << sizes.size() << " above 5x mean, "
              << nCapped << " capped at " << sizeCap
              << " (feedback-divergent)" << std::endl;
    std::ofstream fs("eh_sizes.txt");
    for (const auto s : sizes) fs << s << "\n";
    std::cout << "avalanche sizes written to eh_sizes.txt" << std::endl;
  }

  if (!doMIP) return 0;
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
  unsigned long ncl = 0;
  double yLast = 0.;
  while (track.GetCluster(xc, yc, zc, tc, nc, ec, extra)) {
    if (++ncl % 10 == 0) std::cout << "  cluster " << ncl << std::endl;
    yLast = yc;
    nPrimary += nc;
    for (int k = 0; k < nc; ++k) {
      aval.AvalancheElectronHole(xc, yc, zc, tc);
      std::size_t ne = 0, ni = 0;
      aval.GetAvalancheSize(ne, ni);
      nTotal += ne;
    }
  }
  std::cout << "MIP: " << ncl << " clusters, last at y = "
            << yLast * 1.e4 << " um; " << nPrimary
            << " primary e-h pairs, " << nTotal
            << " electrons after multiplication -> collected-charge gain "
            << double(nTotal) / double(nPrimary) << std::endl;

  double qInt = 0., iMax = 0.;
  for (unsigned int i = 0; i < 800; ++i) {
    const double s = sensor.GetSignal("pad", i);
    qInt += s;
    if (std::abs(s) > std::abs(iMax)) iMax = s;
  }
  std::cout << "signal check: peak bin = " << iMax
            << ", integral (arb.) = " << qInt << std::endl;

  TCanvas c("c", "", 800, 600);
  ViewSignal vs;
  vs.SetSensor(&sensor);
  vs.SetCanvas(&c);
  vs.PlotSignal("pad");
  c.SaveAs("signal_pad.pdf");
  app.Run(true);
}