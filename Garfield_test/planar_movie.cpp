#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <TApplication.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TLatex.h>
#include <TProfile.h>
#include <TROOT.h>
#include <TSystem.h>

#include "Garfield/AvalancheMC.hh"
#include "Garfield/ComponentConstant.hh"
#include "Garfield/ComponentTcad2d.hh"
#include "Garfield/GeometrySimple.hh"
#include "Garfield/MediumSilicon.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/SolidBox.hh"
#include "Garfield/TrackHeed.hh"
#include "Garfield/ViewSignal.hh"
#include "Garfield/ViewDrift.hh"

using namespace Garfield;

int main(int argc, char* argv[]) {
  const std::string file = argc > 1 ? argv[1]
      : "/home/ahaines561/HEP/MAS/Silvaco_dat/lgad.sta";
  const double xTrackUm = argc > 2 ? std::atof(argv[2]) : 20.;
  int rootArgc = 1;
  char* rootArgv[] = {argv[0], nullptr};
  TApplication app("app", &rootArgc, rootArgv);
  gROOT->SetBatch(kTRUE);
  // all plots go into tutorial_plots/ inside the working directory
  const std::string outDir = "tutorial_plots";
  gSystem->mkdir(outDir.c_str(), kTRUE);

  MediumSilicon si;
  si.SetTemperature(293.15);
  si.SetImpactIonisationModelVanOverstraetenDeMan();

  ComponentTcad2d cmp;
  if (!cmp.InitialiseSilvaco(file)) return 1;
  cmp.SetMedium("3", &si);
  cmp.SetRangeZ(-5.e-4, 5.e-4);

  double bx0, by0, bz0, bx1, by1, bz1;
  cmp.GetBoundingBox(bx0, by0, bz0, bx1, by1, bz1);
  const double x0 = xTrackUm * 1.e-4 + 0.13e-4;

  // field profile
  std::vector<double> py, pe;
  double yTop = 1., yBot = -1., eMax = 0., yGain = 0.;
  for (int i = 0; i <= 400; ++i) {
    const double y = by0 + (by1 - by0) * i / 400.;
    double ex, ey, ez, v; int st; Medium* m = nullptr;
    cmp.ElectricField(x0, y, 0., ex, ey, ez, v, m, st);
    if (st != 0) continue;
    if (yTop > yBot) yTop = y;
    yBot = y;
    const double e = std::sqrt(ex * ex + ey * ey);
    py.push_back(y * 1.e4);
    pe.push_back(e * 1.e-3);       // kV/cm
    if (e > eMax) { eMax = e; yGain = y; }
  }
  const double d = yBot - yTop;
  std::cout << "silicon: y = [" << yTop * 1.e4 << ", " << yBot * 1.e4
            << "] um, peak field " << eMax * 1.e-3 << " kV/cm at y = "
            << yGain * 1.e4 << " um" << std::endl;

  TCanvas c1("c1", "", 700, 500);
  TGraph gProf(py.size(), py.data(), pe.data());
  gProf.SetTitle(";depth y [#mum];|E| [kV/cm]");
  gProf.SetLineWidth(2);
  gProf.Draw("AL");
  c1.Print((outDir + "/field_profile.pdf").c_str());

  // weighting field/potential
  ComponentConstant wcmp;
  wcmp.SetArea(bx0, yTop, -5.e-4, bx1, yBot, 5.e-4);
  wcmp.SetMedium(&si);
  wcmp.SetElectricField(0., 0., 0.);
  wcmp.SetWeightingField(0., 1. / d, 0., "pad");
  // Signals are computed from the weighting POTENTIAL: anchor psi = 1
  // at the readout plane; it extends linearly to 0 at the backside.
  wcmp.SetWeightingPotential(0.5 * (bx0 + bx1), yTop, 0., 1.);

  Sensor sensor;
  sensor.AddComponent(&cmp);
  sensor.AddElectrode(&wcmp, "pad");
  sensor.SetTimeWindow(0., 0.005, 800);
  sensor.SetArea(bx0, yTop + 0.02e-4, -5.e-4, bx1, yBot - 0.5e-4, 5.e-4);

  AvalancheMC aval;
  aval.SetSensor(&sensor);
  aval.SetTimeSteps(2.e-3); // 2 ps
  aval.EnableSignalCalculation();
  aval.EnableAvalancheSizeLimit(20000); //breakdown limit

  const double yInj = std::min(yGain + 5.e-4, 0.5 * (yGain + yBot));
  TH1F hAval("hAval", ";avalanche size [electrons];injections", 60, 0.5, 60.5);
  double gE = 0., gEH = 0.;
  const int nInj = 100;
  for (int i = 0; i < nInj; ++i) {
    std::size_t ne = 0, ni = 0;
    aval.AvalancheElectron(x0, yInj, 0., 0.);
    aval.GetAvalancheSize(ne, ni);
    gE += double(ne) / nInj;
    aval.AvalancheElectronHole(x0, yInj, 0., 0.);
    aval.GetAvalancheSize(ne, ni);
    gEH += double(ne) / nInj;
    hAval.Fill(ne);
  }
  std::cout << "G_e = " << gE << ", G_eh = " << gEH << std::endl;

  TCanvas c2("c2", "", 700, 500);
  c2.SetLogy();
  hAval.SetStats(false);
  hAval.SetMinimum(0.5);
  hAval.SetFillColor(kAzure - 9);
  hAval.Draw("hist");
  TLatex tl;
  tl.SetNDC();
  tl.DrawLatex(0.55, 0.80, Form("G_{e}  = %.2f", gE));
  tl.DrawLatex(0.55, 0.72, Form("G_{eh} = %.2f", gEH));
  c2.Print((outDir + "/aval_sizes.pdf").c_str());

  // MIP, charge deposition, gain vs depth, signal
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

  ViewDrift driftView;
  driftView.SetArea(bx0, yTop, -5.e-4, bx1, yBot, 5.e-4);
  track.EnablePlotting(&driftView);

  TH1F hDep("hDep", ";depth y [#mum];primary electrons / #mum",
            50, yTop * 1.e4, yBot * 1.e4);
  TProfile pGain("pGain", ";depth y [#mum];#LTn_{e} after / n_{e} deposited#GT",
                 25, yTop * 1.e4, yBot * 1.e4);
  sensor.ClearSignal();
  track.NewTrack(x0, yTop + 0.01e-4, 0., 0., 0., 1., 0.);
  double xc, yc, zc, tc, ec, extra;
  int nc = 0;
  struct Cl { double x, y, z, t; int n; };
  std::vector<Cl> clusters;                 // kept for the movie
  while (track.GetCluster(xc, yc, zc, tc, nc, ec, extra)) {
    clusters.push_back({xc, yc, zc, tc, nc});
    hDep.Fill(yc * 1.e4, nc / ((yBot - yTop) * 1.e4 / 50.));
    aval.AvalancheElectronHole(xc, yc, zc, tc);
    std::size_t ne = 0, ni = 0;
    aval.GetAvalancheSize(ne, ni);
    pGain.Fill(yc * 1.e4, double(ne) / nc);
  }

  TCanvas c3("c3", "", 700, 500);
  hDep.SetFillColor(kOrange - 9);
  hDep.Draw("hist");
  c3.Print((outDir + "/charge_deposit.pdf").c_str());

  TCanvas c4("c4", "", 700, 500);
  pGain.SetMinimum(0.);
  pGain.SetMarkerStyle(20);
  pGain.Draw();
  c4.Print((outDir + "/gain_vs_depth.pdf").c_str());

  // animated drift + live signal (pattern: Examples/Silicon/planar_movie.C)
  // The MIP's carriers are re-queued and transported in resumable time
  // slices; each slice appends a GIF frame with the drift lines so far
  // (left) and the signal accumulated so far (right).
  Sensor movieSensor;
  movieSensor.AddComponent(&cmp);
  movieSensor.AddElectrode(&wcmp, "pad");
  movieSensor.SetTimeWindow(0., 0.005, 800);
  // deliberately NO SetArea fence: cluster carriers sit at the silicon
  // edges and must be allowed to drift.
  movieSensor.ClearSignal();
  AvalancheMC drift(&movieSensor);
  drift.SetDistanceSteps(1.e-4);            // 1 um steps: smooth lines
  drift.EnableSignalCalculation();
  drift.EnablePlotting(&driftView);
  for (const auto& cl : clusters) {
    for (int j = 0; j < cl.n; ++j) {
      drift.AddElectron(cl.x, cl.y, cl.z, cl.t);
      drift.AddHole(cl.x, cl.y, cl.z, cl.t);
    }
  }

  ViewSignal signalView;
  signalView.SetSensor(&movieSensor);

  TCanvas c5("c5", "", 1400, 600);
  c5.Divide(2, 1);
  driftView.SetCanvas((TPad*)c5.cd(1));
  signalView.SetCanvas((TPad*)c5.cd(2));

  double t0 = 0.;
  const double dt = 0.02;                   // ns per frame
  const std::size_t nFrames = 80;           // 1.6 ns total
  for (std::size_t f = 0; f < nFrames; ++f) {
    driftView.Clear();
    drift.SetTimeWindow(t0, t0 + dt);
    drift.ResumeAvalanche();
    c5.cd(1);
    driftView.Plot2d(true, true);
    c5.cd(2);
    signalView.PlotSignal("pad");
    if (f == nFrames - 1) {
      c5.Print((outDir + "/drift_anim.gif++").c_str());
    } else {
      c5.Print((outDir + "/drift_anim.gif+3").c_str());
    }
    t0 += dt;
  }
  drift.UnsetTimeWindow();

  // Static plots from the finished movie state: all drift lines, and
  // the complete signal.
  TCanvas c6("c6", "", 600, 600);
  driftView.SetCanvas(&c6);
  driftView.Plot2d(true, false);
  c6.Print((outDir + "/drift_lines.pdf").c_str());

  TCanvas c7("c7", "", 700, 500);
  signalView.SetCanvas(&c7);
  signalView.PlotSignal("pad");
  c7.Print((outDir + "/signal_pad.pdf").c_str());

  std::cout << "wrote " << outDir << "/: field_profile.pdf aval_sizes.pdf "
            << "charge_deposit.pdf gain_vs_depth.pdf drift_lines.pdf "
            << "drift_anim.gif signal_pad.pdf" << std::endl;
  return 0;
}