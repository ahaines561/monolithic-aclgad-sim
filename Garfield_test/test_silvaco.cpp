#include <iostream>
#include <iomanip>
#include <TApplication.h>
#include <TCanvas.h>

#include "Garfield/ComponentTcad2d.hh"
#include "Garfield/MediumSilicon.hh"
#include "Garfield/ViewField.hh"

using namespace Garfield;
int main(int argc, char* argv[]) {
  const std::string file =
      argc > 1 ? argv[1]
               : "/home/ahaines561/HEP/MAS/Silvaco_dat/"
                 "30kev5umsp_5um_X26.sta";

  ComponentTcad2d cmp;
  cmp.EnableDebugging();
  if (!cmp.InitialiseSilvaco(file)) {
    std::cerr << "InitialiseSilvaco FAILED on " << file << "\n";
    return 1;
  }
  cmp.DisableDebugging();

  MediumSilicon si;
  for (std::size_t i = 0; i < cmp.GetNumberOfRegions(); ++i) {
    cmp.SetMedium(i, &si);
  }
  // Bounding box check
  double xmin, ymin, zmin, xmax, ymax, zmax;
  if (cmp.GetBoundingBox(xmin, ymin, zmin, xmax, ymax, zmax)) {
    std::cout << "Bounding box [cm]: x[" << xmin << ", " << xmax
              << "]  y[" << ymin << ", " << ymax << "]\n";
  }
 // for mesh visualization
  std::cout << std::scientific << std::setprecision(4);
  std::cout << "\n  y[um]      Ex[V/cm]      Ey[V/cm]   status region\n";
  const double x = 125.e-4;              // 125 um in cm
  for (double yum = 49.5; yum > 0.; yum -= 2.0) {
    const double y = yum * 1.e-4;
    double ex, ey, ez, v;
    Medium* m = nullptr;
    int status = 0;
    cmp.ElectricField(x, y, 0., ex, ey, ez, v, m, status);
    std::cout << std::fixed << std::setprecision(2) << std::setw(7) << yum
              << std::scientific << std::setprecision(4)
              << "  " << std::setw(12) << ex
              << "  " << std::setw(12) << ey
              << std::fixed << "   " << status << "\n";
  }
  auto probe = [&](double yum) {
    double ex, ey, ez, v; Medium* m = nullptr; int status = 0;
    cmp.ElectricField(125.e-4, yum * 1.e-4, 0., ex, ey, ez, v, m, status);
    std::cout << "  y=" << yum << " um : Ey=" << ey << " V/cm  V=" << v
              << " status=" << status << "\n";
  };
  std::cout << "\nspot checks:\n";
  std::cout << std::scientific << std::setprecision(4);
  probe(25.0);
  probe(1.0);

  TApplication app("app", &argc, argv);

TCanvas* c1 = new TCanvas("c1", "Silvaco Ey", 900, 500);
ViewField fview;
fview.SetComponent(&cmp);
fview.SetCanvas(c1);
fview.SetArea(0., 0., 250.e-4, 50.e-4);
fview.SetNumberOfSamples2d(400, 200);
fview.Plot("ey", "colz"); 
c1->Update();

  return 0;
}