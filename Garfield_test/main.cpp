#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <TApplication.h>
#include <TCanvas.h>

#include "Garfield/MediumSilicon.hh"
#include "Garfield/ComponentTcad2d.hh"
#include "Garfield/ComponentUser.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/TrackHeed.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/ViewSignal.hh"

using namespace Garfield;

int main(int argc, char* argv[]) {
    TApplication app("app", &argc, argv);
    //open file
  const std::string file =
      argc > 1 ? argv[1]
               : "/home/ahaines561/HEP/MAS/Silvaco_dat/"
                 "30kev5umsp_5um_X26.sta";

    MediumSilicon si;
    si.SetTemperature(293.15);
    si.SetImpactIonisationModelOkutoCrowell();

    ComponentTcad2d field;
    if (!field.InitialiseSilvaco(file)) {
        std::cerr << "InitialiseSilvaco FAILED\n";
        return 1;
    }
    // Bind the Silicon medium using the direct integer Index from PrintRegions
    for (int i = 1; i <= 13; ++i) {
        field.SetMedium(i, &si);
    }
    field.PrintRegions();

    constexpr double yBot = -3.60e-4;
    constexpr double ySi_Top = 0.0e-4;
    constexpr double ySi_Bot = 50.50e-4;

    ComponentUser wpad;
    // SetArea: (xMin, yMin, zMin, xMax, yMax, zMax)
    wpad.SetArea(0., yBot, -1., 250.e-4, ySi_Bot, 1.);
    wpad.SetMedium(&si);
    
    wpad.SetWeightingPotential(
        [ySi_Top, ySi_Bot](const double, const double y, const double){
            if (y < ySi_Top) return 1.0;
            if (y > ySi_Bot) return 0.0;
            return (ySi_Bot - y) / (ySi_Bot - ySi_Top);
        }, "pad"); 
    wpad.SetWeightingField(
        [ySi_Top, ySi_Bot](const double, const double y, const double, double& wx, double& wy, double& wz) {
            wx = 0.; 
            wy = (y >= ySi_Top && y <= ySi_Bot) ? 1.0 / (ySi_Bot - ySi_Top) : 0.; 
            wz = 0.;
        }, "pad");

    // Sensor
    Sensor sensor;
    sensor.AddComponent(&field);
    sensor.AddElectrode(&wpad, "pad");
    
    // Set a 3 ns time window with 1000 bins (3 ps resolution)
    const int nBins = 1000;
    const double tmax = 3.0;
    const double dt = tmax / nBins;
    sensor.SetTimeWindow(0., dt, nBins);
    sensor.SetArea(0., yBot, -1., 250.e-4, ySi_Bot, 1.);

    // MIP
    TrackHeed track(&sensor);
    track.SetParticle("electron");
    track.SetMomentum(9.e10);
    track.EnableDeltaElectronTransport();

    // Charge Transport and Avalanche
    AvalancheMC aval(&sensor);
    aval.SetDistanceSteps(0.1e-4);
    aval.EnableSignalCalculation();
    aval.EnableAvalancheSizeLimit(100000);
    aval.SetTimeWindow(0., 10.);

    // run the Simulation
    std::cout << "Starting MIP Simulation \n";
    std::cout << ". \n";
    std::cout << ". \n";
    std::cout << ". \n";
    
    // Shoot the track straight up the middle of the 250um wide sensor
    track.NewTrack(118.e-4, 45.e-4, 0., 0., 0., -1., 0.);

    // capture the primary ionization once so we can cap how many we transport
    struct Seed { double x, y, z, t; };
    std::vector<Seed> seeds;
    for (const auto& cluster : track.GetClusters()) {
        for (const auto& e : cluster.electrons) {
            seeds.push_back({e.x, e.y, e.z, e.t});
        }
    }
    const std::size_t nUse = std::min<std::size_t>(seeds.size(), 50);

    // Drift the primary electron and calculate multiplication in the gain layer
    int totElectrons = 0, totHoles = 0, eDone = 0;
    for (std::size_t i = 0; i < nUse; ++i) {
        const auto& s = seeds[i];
        aval.AvalancheElectronHole(s.x, s.y, s.z, s.t);
        totElectrons += aval.GetElectrons().size();
        totHoles += aval.GetHoles().size();
        if (++eDone % 10 == 0) {
            std::cout << "  Processed " << eDone << "/" << nUse 
                      << " primary electrons\r" << std::flush;
        }
    }
    std::cout << "\nTrack finished. Primary pairs: " << seeds.size() << "\n";
    std::cout << "Transported: " << nUse << "\n";
    std::cout << "Total electrons after avalanche: " << totElectrons << "\n";
    std::cout << "Total holes after avalanche: " << totHoles << "\n";
    if (nUse > 0) {
        std::cout << "Carrier-count gain: "
                  << static_cast<double>(totElectrons) / nUse << "\n";
    }

    // Plotting
    double qPad = 0.0;
    for (int i = 0; i < nBins; ++i) {
        qPad += sensor.GetSignal("pad", i) * dt;
    }
    std::cout << "Integrated charge on pad: " << qPad << " fC\n";

    TCanvas* c1 = new TCanvas("c1", "AC-LGAD Signal", 800, 600);
    ViewSignal sv(&sensor);
    sv.SetCanvas(c1);
    sv.PlotSignal("pad", "t");
    c1->Update();
    c1->SaveAs("silvaco_signal.png");

    app.Run();
    return 0;
}