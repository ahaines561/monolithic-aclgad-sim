#include <iostream>
#include <vector>

#include <TApplication.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TH2F.h>

#include "Garfield/ComponentTcad2d.hh"
#include "Garfield/MediumSilicon.hh"

using namespace Garfield;

int main(int argc, char* argv[]) {
  TApplication app("app", &argc, argv);

  const std::string file =
      argc > 1 ? argv[1]
               : "/home/ahaines561/HEP/MAS/Silvaco_dat"
                 "30kev5umsp_5um_X26.sta";

  ComponentTcad2d cmp;
  if (!cmp.InitialiseSilvaco(file)) {
    std::cerr << "InitialiseSilvaco FAILED\n";
    return 1;
  }

  TCanvas* c = new TCanvas("c", "Silvaco mesh", 1200, 350);
  TH2F* frame = new TH2F("frame", "Silvaco mesh;x [um];y [um]",
                         10, 0., 250., 10, -4., 51.);
  frame->SetStats(0);
  frame->Draw();

  const std::size_t nEl = cmp.GetNumberOfElements();
  std::vector<TLine*> lines;
  lines.reserve(nEl * 3);

  for (std::size_t i = 0; i < nEl; ++i) {
    std::vector<std::size_t> nodes;
    if (!cmp.GetElementNodes(i, nodes, false)) continue;
    if (nodes.size() < 3) continue;
    double x[3], y[3], z;
    bool ok = true;
    for (int k = 0; k < 3; ++k) {
      if (!cmp.GetNode(nodes[k], x[k], y[k], z)) { ok = false; break; }
      x[k] *= 1.e4;   // cm -> um
      y[k] *= 1.e4;
    }
    if (!ok) continue;
    for (int k = 0; k < 3; ++k) {
      const int k2 = (k + 1) % 3;
      TLine* L = new TLine(x[k], y[k], x[k2], y[k2]);
      L->SetLineColor(kAzure + 1);
      L->SetLineWidth(1);
      L->Draw();
      lines.push_back(L);
    }
  }

  c->Update();
  c->SaveAs("silvaco_mesh.png");
  std::cout << "Drew " << nEl << " elements. Saved silvaco_mesh.png\n";
  std::cout << "Close the canvas or Ctrl-C to exit.\n";
  app.Run();
  return 0;
}