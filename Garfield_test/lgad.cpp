#include <cmath>
#include <iostream>
#include <TApplication.h>
#include <TCanvas.h>
#include "Garfield/ComponentTcad2d.hh"
#include "Garfield/ComponentConstant.hh"
#include "Garfield/MediumSilicon.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/TrackHeed.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/ViewSignal.hh"

using namespace Garfield;

int main(int argc, char* argv[]) {
  TApplication app("app", &argc, argv);
  const std::string file = argc > 1 ? argv[1]
      : "/home/ahaines561/HEP/MAS/Silvaco_dat/30kev5umsp_5um_X26.sta";

  // Device extent [cm]: 250 x 50.5 um, silicon from y = 0 to 50 um.
  constexpr double xDev = 250.e-4;
  constexpr double yTop = 0.;
  constexpr double yBot = 50.e-4;
  constexpr double d = yBot - yTop;      // drift thickness
  const double x0 = 0.5 * xDev;          // probe/track x (strip centre?)

  ComponentTcad2d cmp;
  if (!cmp.InitialiseSilvaco(file)) return 1;

  MediumSilicon si;
  si.SetTemperature(293.15);
  si.SetImpactIonisationModelVanOverstraetenDeMan();

  cmp.SetMedium(0, &si);
  for (std::size_t i = 1; i < 13; ++i) cmp.UnsetDriftRegion(i);
  cmp.SetRangeZ(-5.e-4, 5.e-4);

  std::cout << "# y[um]   V[V]   Ey[V/cm]\n";
  double eMax = 0., yGain = 0.;
  for (int i = 0; i <= 200; ++i) {
    const double y = yTop + (d * i) / 200;
    double ex, ey, ez, v; int st; Medium* m = nullptr;
    cmp.ElectricField(x0, y, 0., ex, ey, ez, v, m, st);
    if (st != 0) continue;
    const double e = std::sqrt(ex * ex + ey * ey);
    if (e > eMax) { eMax = e; yGain = y; }
    std::cout << y * 1.e4 << "  " << v << "  " << ey << "\n";
  }
  std::cout << "Peak field " << eMax << " V/cm at y = " << yGain * 1.e4
            << " um  (gain layer)\n";

  ComponentConstant wcmp;
  wcmp.SetArea(0., yTop, -5.e-4, xDev, yBot, 5.e-4);
  wcmp.SetMedium(&si);
  wcmp.SetElectricField(0., 0., 0.);
  wcmp.SetWeightingField(0., 1. / d, 0., "pad");

  Sensor sensor;
  sensor.AddComponent(&cmp);
  sensor.AddElectrode(&wcmp, "pad");
  sensor.SetTimeWindow(0., 0.005, 800);

  AvalancheMC aval;
  aval.SetSensor(&sensor);
  aval.EnableSignalCalculation();
  aval.SetTimeSteps(5.e-4);

  const int nInject = 50;
  double sumGain = 0.;
  for (int i = 0; i < nInject; ++i) {
    if (i % 10 == 0) std::cout << "  injection " << i << std::endl;
    const double xi = x0 + (i % 5 - 2) * 2.e-4;
    aval.AvalancheElectron(xi, yGain + 5.e-4, 0., 0.);
    std::size_t ne = 0, ni = 0;
    aval.GetAvalancheSize(ne, ni);
    sumGain += ne;
  }
  std::cout << "Single-electron gain (mean of " << nInject
            << " injections): " << sumGain / nInject << std::endl;

  TrackHeed track;
  track.SetSensor(&sensor);
  track.SetParticle("pi");
  track.SetMomentum(180.e9);
  sensor.ClearSignal();

  track.NewTrack(x0, yTop + 0.01e-4, 0., 0., 0., 1., 0.);
  double xc, yc, zc, tc, ec, extra;
  int nc = 0;
  unsigned long nPrimary = 0, nTotal = 0;
  unsigned long ncl = 0;
  while (track.GetCluster(xc, yc, zc, tc, nc, ec, extra)) {
    if (++ncl % 200 == 0) std::cout << "  cluster " << ncl << std::endl;
    nPrimary += nc;
    for (int k = 0; k < nc; ++k) {
      aval.AvalancheElectron(xc, yc, zc, tc);
      std::size_t ne = 0, ni = 0;
      aval.GetAvalancheSize(ne, ni);
      nTotal += ne;
      aval.DriftHole(xc, yc, zc, tc);
    }
  }
  std::cout << "MIP: " << nPrimary << " primary e-h pairs, " << nTotal
            << " electrons after multiplication -> collected-charge gain "
            << double(nTotal) / double(nPrimary) << std::endl;

  TCanvas c("c", "", 800, 600);
  ViewSignal vs;
  vs.SetSensor(&sensor);
  vs.SetCanvas(&c);
  vs.PlotSignal("pad");
  c.SaveAs("signal_pad.pdf");
  app.Run(true);
}