#include <TApplication.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "TGraph2D.h"

#include "Garfield/MediumSilicon.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/SolidBox.hh"
#include "Garfield/GeometrySimple.hh"
#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/TrackHeed.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/ViewGeometry.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/ViewDrift.hh"
#include "Garfield/ViewMedium.hh" 
#include "Garfield/ComponentUser.hh"
#include "Garfield/ViewSignal.hh"
#include "Garfield/ComponentTcad3d.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/ComponentConstant.hh"
#include "Garfield/ComponentTcad2d.hh"
#include "Garfield/ComponentGrid.hh"
#include "Garfield/DriftLineRKF.hh"
#include "Garfield/ViewDrift.hh"

using namespace Garfield;

//helper function to create interpolation map from unstructured silvaco tcad data
TGraph2D* LoadField(const std::string& grd, const std::string& dat, int col) {
    auto* g = new TGraph2D();
    std::ifstream fg(grd), fd(dat);
    if (!fg.is_open()) { std::cerr << "Cannot open " << grd << "\n"; return g; }
    if (!fd.is_open()) { std::cerr << "Cannot open " << dat << "\n"; return g; }
    double xu, yu, v1, v2;
    int n = 0;
    while (fg >> xu >> yu) {
        fd >> v1 >> v2;                       // .dat always has 2 cols (Ex, Ey)
        const double xc = xu * 1.e-4;         // um -> cm
        const double yc = yu * 1.e-4;
        g->SetPoint(n, xc, yc, (col == 1) ? v1 : v2);
        ++n;
    }
    std::cout << "[LoadField] " << n << " nodes from " << grd << "\n";
    return g;
}

int main (int argc, char *argv[]) {
//grid geometry from silvaco tcad data
constexpr double xMin = 0.;
constexpr double xMax = 250.e-4;
constexpr double yMin = -3.6e-4; 
constexpr double yMax = 50.5e-4;
constexpr double zHalf = 1.e-4;
constexpr double dDepth = yMax - yMin;
//coordinates
double x0 = 125.e-4, y0 = 50.0e-4, z0 = 0., t0 = 0.;
double dx = 0., dy = -1., dz = 0.;
int tot_Electrons = 0;
int tot_Holes = 0;

    const std::string dir = "/home/ahaines561/HEP/Si_Test/data_folder/";
    TGraph2D* gEx = LoadField(dir + "k30_ElectricField.grd", dir + "k30_ElectricField.dat", 1);
    TGraph2D* gEy = LoadField(dir + "k30_ElectricField.grd", dir + "k30_ElectricField.dat", 2);

if (gEx->GetN() == 0 || gEy->GetN() == 0) {
    std::cerr << "ERROR: no field points loaded — check dir path and filenames\n";
    return 1;
}
else{
    std::cout << "eGx and eGy loaded fine" << std::endl;
}

TApplication app("app", &argc, argv);
MediumSilicon si;
si.SetTemperature(293.15);
si.SetImpactIonisationModelOkutoCrowell();
//p-type doping
si.SetDopingMobilityModelMasetti();
si.SetDoping('p', 1.0e12);

//load Silvaco data
ComponentUser field;
field.SetArea(xMin, yMin, -zHalf, xMax, yMax, zHalf);
field.SetMedium(&si);
field.SetElectricField(
    [gEx, gEy](const double x, const double y, const double /*z*/,
               double& ex, double& ey, double& ez) {
        ex = gEx->Interpolate(x, y);
        ey = gEy->Interpolate(x, y);
        ez = 0.;
    });

// tcad.SetWeightingField("aclgad_WeightingPotential.grd", "aclgad_WeightingPotential.dat", "pad");
// sensor.AddElectrode(&tcad, "pad"); //weighted field is needed

ComponentUser wpad;
wpad.SetArea(xMin, yMin, -zHalf, xMax, yMax, zHalf);
wpad.SetMedium(&si);
//weighting potential and field for the pad electrode
wpad.SetWeightingPotential(
    [yMin, yMax](const double, const double y, const double) {
        return (yMax - y) / (yMax - yMin);   // W=1 at bottom pad, W=0 at top
    }, "pad");
wpad.SetWeightingField(
    [yMax, yMin](const double, const double, const double,
                 double& wx, double& wy, double& wz) {
        wx = 0.; wy = -1. / (yMax - yMin); wz = 0.;   // note sign
    }, "pad");

Sensor sensor;
sensor.AddComponent(&field);
sensor.AddElectrode(&wpad, "pad");
const std::size_t nTimeBins = 1000;
const double tmin = 0.;
const double tmax = 5.; 
const double tstep = (tmax - tmin) / nTimeBins;
sensor.SetTimeWindow(tmin, tstep, nTimeBins);

field.SetElectricField(
    [gEx, gEy](const double x, const double y, const double,
               double& ex, double& ey, double& ez) {
        ex = 0.;                      // TEMP test
        ey = gEy->Interpolate(x, y);
        ez = 0.;
    });

for (double yy = yMax; yy > 0; yy -= 2.e-4) {
    double ex, ey, ez; Medium* m; int status;
    sensor.ElectricField(125.e-4, yy, 0., ex, ey, ez, m, status);
    std::cout << "y=" << yy*1e4 << "um  Ey=" << ey
              << "  status=" << status << "\n";
}

//MIP
TrackHeed track(&sensor);
track.SetParticle("pion");
track.SetMomentum(180.e9);
//enforce delta ray degradation in si
track.EnableDeltaElectronTransport();

//AvalancheMC 
AvalancheMC aval(&sensor);
aval.SetDistanceSteps(1.e-4);   // 1 μm per step [cm]
aval.EnableSignalCalculation();

//runge-kutta-fehlberg integration
//drift
ViewDrift driftView;
DriftLineRKF drift(&sensor);
drift.EnableSignalCalculation();
drift.EnablePlotting(&driftView);
driftView.SetArea(xMin, yMin, -zHalf, xMax, yMax, zHalf);
track.NewTrack(x0, y0, z0, t0, dx ,dy, dz);
track.EnablePlotting(&driftView);
aval.EnablePlotting(&driftView);

//the avalance simulation
for (const auto& cluster : track.GetClusters()) {
    for (const auto& e : cluster.electrons) {
        //Electron drift + LGAD avalanche gain
        aval.AvalancheElectron(e.x, e.y, e.z, e.t, 0.1);
        tot_Electrons += static_cast<int>(aval.GetElectrons().size());
        drift.DriftHole(e.x, e.y, e.z, e.t); ++tot_Holes;
        }
    }

//visualize
TCanvas* c1 = new TCanvas("c1", "lgad signal", 800, 600);
ViewSignal signalView(&sensor);
signalView.SetCanvas(c1);
signalView.PlotSignal("pad", "t");
c1->Update();

TCanvas* c2 = new TCanvas("c2", "LGAD Drift Lines", 800, 600);
driftView.SetCanvas(c2);
driftView.Plot(true, false);
c2->Update();

app.Run();
return 0;
}