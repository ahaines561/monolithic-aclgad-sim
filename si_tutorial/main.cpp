#include <TApplication.h>
#include <TCanvas.h>
#include <TGraph2D.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include "Garfield/MediumSilicon.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/ComponentUser.hh"
#include "Garfield/TrackHeed.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/ViewDrift.hh"
#include "Garfield/ViewSignal.hh"
#include "Garfield/ViewMedium.hh"
#include "Garfield/ViewField.hh"

using namespace Garfield;

//load field components from Silvaco .grd + .dat 
TGraph2D* LoadField(const std::string& grd, const std::string& dat, int col) {
    auto* g = new TGraph2D();
    std::ifstream fg(grd), fd(dat);
    if (!fg.is_open()) { std::cerr << "Cannot open " << grd << "\n"; return g; }
    if (!fd.is_open()) { std::cerr << "Cannot open " << dat << "\n"; return g; }
    double xu, yu, v1, v2;
    int n = 0;
    while (fg >> xu >> yu) {
        fd >> v1 >> v2;
        g->SetPoint(n, xu * 1.e-4, yu * 1.e-4, (col == 1) ? v1 : v2);  // um -> cm
        ++n;
    }
    std::cout << "[LoadField] " << n << " nodes from " << grd << "\n";
    return g;
}

int main(int argc, char* argv[]) {
    TApplication app("app", &argc, argv);
        const std::string dir = "/home/ahaines561/HEP/Si_Test/data_folder/";
    constexpr double xMin = 0., xMax = 250.e-4;
    constexpr double yMin = -3.6e-4,yMax = 50.5e-4;
    constexpr double zHalf = 1.e-4;
    constexpr double yBot = 0.0e-4;
    constexpr double yTop = 49.5e-4;
    const int nTracks = 200;
    double qSum = 0.0, qSum2 = 0.0;
    double qPad = 0.0;
    const int nBins = 1000;
    const double tmax = 3.0;
    const double dt = tmax / nBins;
    double lnGain = 0.0;
    const double dy = 0.01e-4;
    int nClusters = 0, nPrimaryE = 0;
    int totElectrons = 0, totHoles = 0;
    int eDone = 0;

    //si medium
    MediumSilicon si;
    si.SetImpactIonisationModelOkutoCrowell();
    si.SetTemperature(300.);
    si.SetDoping('p', 1.0e12);

    //load tcad file
    TGraph2D* gEx = LoadField(dir + "k30_ElectricField.grd",
                              dir + "k30_ElectricField.dat", 1);
    TGraph2D* gEy = LoadField(dir + "k30_ElectricField.grd",
                              dir + "k30_ElectricField.dat", 2);
    if (gEy->GetN() == 0) {
        std::cerr << "ERROR: no field points loaded — check path/filenames\n";
        return 1;
    }

    //drift field
    ComponentUser field;
    field.SetArea(xMin, yMin, -zHalf, xMax, yMax, zHalf);
    field.SetMedium(&si);
    field.SetElectricField(
    [gEx, gEy](const double x, const double y, const double,
               double& ex, double& ey, double& ez) {
        ex = gEx->Interpolate(x, y);
        ey = gEy->Interpolate(x, y);
        ez = 0.;
    });

    //weighting potential and field for the pad electrode
    ComponentUser wpad;
    wpad.SetArea(xMin, yMin, -zHalf, xMax, yMax, zHalf);
    wpad.SetMedium(&si);
    wpad.SetWeightingPotential(
        [yBot, yTop](const double, const double y, const double){
            return (yTop - y) / (yTop - yBot);
        }, "pad");
    wpad.SetWeightingField(
        [yBot, yTop](const double, const double, const double,
                     double& wx, double& wy, double& wz) {
            wx = 0.; wy = -1. / (yTop - yBot); wz = 0.;
        }, "pad");

    //sensor
    Sensor sensor;
    sensor.AddComponent(&field);
    sensor.AddElectrode(&wpad, "pad");
    sensor.SetTimeWindow(0., dt, nBins);

    //sanity check scan field along x
    std::cout <<"Field scan along x=125um \n";
    for (double yy = yTop; yy > yBot; yy -= 2.e-4){
        double ex = 0., ey = 0., ez = 0.;
        Medium* m = nullptr; int status = -99;
        sensor.ElectricField(125.e-4, yy, 0., ex, ey, ez, m, status);
        std::cout << "y=" << yy*1e4 << " um  Ex=" << ex
                  << "  Ey=" << ey << "  status=" << status << "\n";
    }

    //MIP
    TrackHeed track(&sensor);
    track.SetParticle("pion");
    track.SetMomentum(180.e9);
    track.EnableDeltaElectronTransport();

    //drift, no avalanche
    AvalancheMC aval(&sensor);
    aval.SetDistanceSteps(0.1e-4);
    aval.EnableSignalCalculation();

    ViewDrift driftView;
    driftView.SetArea(xMin, yBot, -zHalf, xMax, yTop, zHalf);
    track.EnablePlotting(&driftView);
    aval.EnablePlotting(&driftView);

    //MIP simulation
    track.NewTrack(125.e-4, yTop, 0., 0.,
                   0., -1., 0.);

    for (const auto& cluster : track.GetClusters()) {
        ++nClusters;;
        nPrimaryE += static_cast<int>(cluster.electrons.size());
    }
    std::cout << "Track produced " << nClusters << " clusters, "
              << nPrimaryE << " primary electrons\n";
    if (nClusters == 0) {
        std::cerr << "WARNING: track made no clusters — check track start "
                     "position / direction / momentum.\n";
    }

    for (const auto& cluster : track.GetClusters()) {
        for (const auto& e : cluster.electrons) {
            aval.DriftElectron(e.x, e.y, e.z, e.t);
            totElectrons += static_cast<int>(aval.GetElectrons().size());
            aval.DriftHole(e.x, e.y, e.z, e.t);
            totHoles += static_cast<int>(aval.GetHoles().size());

            if (++eDone % 50 == 0)
                std::cout << "  processed " << eDone << "/" << nPrimaryE
                          << " primary electrons\r" << std::flush;
        }
    }
    std::cout << "Electron endpoints: " << totElectrons
              << "   Holes drifted: " << totHoles << "\n";

    //collected charge on pad: integrate the induced current over the time bins
    for (int i = 0; i < nBins; ++i) {
        qPad += sensor.GetSignal("pad", i) * dt;
    }
    std::cout << "Integrated charge on pad: " << qPad << " fC\n";

    for (double yy = yTop; yy > yBot; yy -= dy) {
    double ex, ey, ez; Medium* m; int status;
    sensor.ElectricField(125.e-4, yy, 0., ex, ey, ez, m, status);
    const double emag = std::sqrt(ex*ex + ey*ey + ez*ez);
    double alpha = 0.0;
    si.ElectronTownsend(ex, ey, ez, 0, 0, 0, alpha);   // alpha in 1/cm
    lnGain += alpha * dy;
    }
    std::cout << "Analytic gain estimate exp(integral alpha dy) = "
            << std::exp(lnGain) << "\n";

    // for (const auto& cluster : track.GetClusters()) {
    //     for (const auto& e : cluster.electrons) {
    //         //electron avalanche (gain); 'true' = also drift the holes it creates
    //         aval.AvalancheElectron(e.x, e.y, e.z, e.t, true);
    //         totElectrons += static_cast<int>(aval.GetElectrons().size());
    //         totHoles      += static_cast<int>(aval.GetHoles().size());
    //     }
    // }

    TCanvas* c1 = new TCanvas("c1", "LGAD signal", 800, 600);
    ViewSignal sv(&sensor);
    sv.SetCanvas(c1);
    sv.PlotSignal("pad", "t");
    c1->Update();

    TCanvas* c2 = new TCanvas("c2", "LGAD Drift Lines", 800, 600);
    driftView.SetCanvas(c2);
    driftView.Plot(true, true);
    c2->Update();

    TCanvas* c3 = new TCanvas("c3", "Drift velocity", 800, 600);
    ViewMedium medView(&si);
    medView.SetCanvas(c3);
    medView.PlotElectronVelocity('e');
    medView.PlotHoleVelocity('e', true);
    c3->Update();

    TCanvas* c4 = new TCanvas("c4", "Weighting potential", 800, 600);
    ViewField wView;
    wView.SetComponent(&wpad);
    wView.SetCanvas(c4);
    wView.SetArea(xMin, yBot, xMax, yTop);
    wView.PlotContourWeightingField("pad", "v");
    c4->Update();

    TCanvas* c5 = new TCanvas("c5", "Ey map", 800, 600);
    ViewField fView;
    fView.SetComponent(&field);
    fView.SetCanvas(c5);
    fView.SetArea(xMin, yBot, xMax, yTop);
    fView.SetNumberOfSamples2d(200, 200);
    fView.Plot("ey", "COLZ");
    c5->Update();

    app.Run();
    return 0;
}