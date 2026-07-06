#include <iostream>
#include <TApplication.h>
#include <TCanvas.h>
#include <TStyle.h>

#include "Garfield/ComponentTcad2d.hh"
#include "Garfield/MediumSilicon.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/ViewFEMesh.hh"

using namespace Garfield;

int main(int argc, char* argv[]) {
  TApplication app("app", &argc, argv);

  const std::string file =
      argc > 1 ? argv[1]
               : "/home/ahaines561/HEP/MAS/Silvaco_dat/"
                 "lgad.sta";

  ComponentTcad2d cmp;
  cmp.InitialiseSilvaco(file);
  MediumSilicon si;
  for (std::size_t i = 0; i < cmp.GetNumberOfRegions(); ++i) {
    cmp.SetMedium(i, &si);
  }

  // the physical silicon for 30keV
  // const double x0 = 0., y0 = -4e-4;
  // const double x1 = 250.e-4, y1 = 50.e-4;
  // the physical silicon for lgad.sta
  const double x0 = 0., y0 = -1.15308e-4;
  const double x1 = 100.e-4, y1 = 50.5e-4;
  

  // Ey colour map
  TCanvas* c1 = new TCanvas("c1", "Silvaco Ey [V/cm]", 900, 450);
  ViewField fEy;
  fEy.SetComponent(&cmp);
  fEy.SetCanvas(c1);
  fEy.SetArea(x0, y0, x1, y1);
  fEy.SetNumberOfSamples2d(400, 200);
  fEy.Plot("ey", "colz");
  c1->Update();

  // |E| colour map (gain layer stands out)
  TCanvas* c2 = new TCanvas("c2", "Silvaco |E| [V/cm]", 900, 450);
  ViewField fEm;
  fEm.SetComponent(&cmp);
  fEm.SetCanvas(c2);
  fEm.SetArea(x0, y0, x1, y1);
  fEm.SetNumberOfSamples2d(400, 200);
  fEm.Plot("e", "colz");
  c2->Update();

  // electrostatic potential
  TCanvas* c3 = new TCanvas("c3", "Silvaco potential [V]", 900, 450);
  ViewField fV;
  fV.SetComponent(&cmp);
  fV.SetCanvas(c3);
  fV.SetArea(x0, y0, x1, y1);
  fV.SetNumberOfSamples2d(400, 200);
  fV.Plot("v", "colz");
  c3->Update();

  // triangular mesh, coloured by region (DOESNT WORK)
  // TCanvas* c4 = new TCanvas("c4", "Silvaco mesh", 900, 450);
  // ViewFEMesh* mesh = new ViewFEMesh(&cmp);
  // mesh->SetCanvas(c4);
  // mesh->SetArea(x0, y0, -1., x1, y1, 1.);
  // mesh->SetPlane(0, 0, 1, 0, 0, 0);
  // mesh->SetFillMesh(false);
  // mesh->EnableAxes();
  // mesh->Plot();
  // c4->Update();

  std::cout << "\nClose the canvases or Ctrl-C to exit.\n";
  app.Run();
  return 0;
}