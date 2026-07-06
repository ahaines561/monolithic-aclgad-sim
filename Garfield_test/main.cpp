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
    TApplication app("app", &argc, argv);
    const std::string file = argc > 1 ? argv[1]
        : "/home/ahaines561/HEP/MAS/Silvaco_dat/lgad.sta";

    ComponentTcad2d cmp;
    cmp.InitialiseSilvaco(file);
    MediumSilicon si;
    si.SetTemperature(293.15);
    si.SetImpactIonisationModelMassey();
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
}